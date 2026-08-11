#!/usr/bin/env python3
"""Step 1: ENDF files for one FNS foil, converted to Arrow on this machine.

Reads the foil's composition from the benchmark's own input deck, expands it
into natural isotopes, and puts each of those ENDF evaluations through NJOY
(resonance reconstruction plus Doppler broadening to 294 K). Each comes out as
a ``Fe56.arrow/`` directory of Arrow record batches, which is the form yats
reads.

    convert_to_arrow.py                                  # iron, the default
    convert_to_arrow.py --case W
    convert_to_arrow.py --case SS316                     # an alloy, 6 elements
    convert_to_arrow.py --endf-dir /path/to/endf --library endf-b8.1

Point ``--endf-dir`` at any directory of ENDF files. It is searched
recursively and each file is identified by the nuclide in its own header
rather than by its name, so ``n-Fe056.tendl``, ``n-026_Fe_056.endf`` and
``whatever.dat`` are all found the same way.

With no ``--endf-dir`` it falls back to fetching the library named by
``--library``, either from ``--tarball`` (a local copy of ``TENDL-n.tgz``) or
from the TENDL website. ``tendl-2017`` and ``tendl-2025`` can both be fetched;
see ``data_source.TENDL_TARBALLS``. That download is around 3 GB however few
nuclides are wanted, because TENDL publishes its neutron sublibrary as a single
archive.

Everything a run produces is kept under the name of the data it used, so a
TENDL-2017 build and a TENDL-2025 build sit side by side instead of overwriting
each other.

Three things come out of the chosen library: the reaction cross sections, the
isomeric branching that says which state each product is made in, and the
reaction topology that says which reaction on which parent gives which product.
What is left on endf-b8.1 is the decay data (half-lives, decay modes, decay
energies), because no neutron-reaction library carries it.

Both chain subsections need that decay data anyway, branching to map a
product's nuclear level to its metastable state and the topology to know which
nuclides exist, so the first run fetches the ENDF/B-8.1 decay sublibrary.
``--decay-dir`` uses a copy you already have. ``--no-branching`` skips the
branching, at the cost of being wrong by up to 6x on foils whose decay heat
comes from an isomer, and ``--no-reactions`` leaves the topology on endf-b8.1.
"""

import argparse
import fnmatch
import os
import pathlib
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
import urllib.request

import data_source
import fns_case
import yats

HERE = pathlib.Path(__file__).resolve().parent

# Both chain subsections need decay data: branching to say which nuclear level
# of an (n,2n) product is which metastable state, and the reaction topology to
# know which nuclides exist at all. The sublibrary is a 10.6 MB download, so it
# is fetched rather than made a prerequisite.
DECAY_URL = ("https://www.nndc.bnl.gov/endf-releases/releases/B-VIII.1/"
             "decay/decay-version.VIII.1.tar.gz")
DECAY_LIBRARY = "endf-b8.1"
DECAY_GLOB = "dec-*.endf"
DECAY_ISOMER_GLOB = "dec-*m[0-9].endf"

# ENDF is a fixed-column format: six 11-character fields, then MAT, MF, MT.
_FIELD = [slice(0, 11), slice(11, 22), slice(22, 33), slice(33, 44)]
_MF = slice(70, 72)
_MT = slice(72, 75)
# ZSYMAM, "26-Fe- 56", is the fifth record of MF=1 MT=451: four numeric records
# come first, then the text section opens with it.
_ZSYMAM = 4
_ZSYMAM_RE = re.compile(r"\s*(\d+)\s*-\s*([A-Za-z]{1,2})\s*-\s*(\d+)")


def isotopes_for(elements):
    """The natural isotopes of each element, as {nuclide: abundance percent}.

    Taken from yats rather than a table kept here, so the isotopes converted
    are exactly the ones yats will expand the material into.
    """
    by_element = yats.data.element_nuclides()
    abundance = yats.data.natural_abundance()
    wanted = {}
    for element in elements:
        if element not in by_element:
            raise SystemExit(f"{element} has no natural isotopes in yats' abundance table")
        for nuclide in by_element[element]:
            wanted[nuclide] = abundance.get(nuclide, 0.0) * 100.0
    return wanted


def endf_float(field):
    """Parse an ENDF number, which writes 2.605600+4 for 2.6056e4."""
    field = field.strip()
    if not field:
        return 0.0
    if "e" in field or "E" in field:
        return float(field)
    return float(field[0] + field[1:].replace("+", "e+").replace("-", "e-"))


def identify(path, lines=8):
    """The nuclide an ENDF file holds, like "Fe56", or None if it is not one.

    Reads the MF=1 MT=451 descriptive section that every ENDF evaluation opens
    with. The name comes from its ZSYMAM field rather than from ZA, so no
    periodic table is needed here, and the isomeric state comes from LIS0 in
    the second record. Metastable evaluations return None: they are separate
    evaluations of the same ZA and are not part of a natural composition. Only
    the first few lines are read, so scanning thousands of files stays cheap.
    """
    try:
        with open(path, "r", errors="replace") as handle:
            head = [next(handle, "") for _ in range(lines)]
    except OSError:
        return None

    for index, line in enumerate(head):
        if len(line) < 75 or line[_MF].strip() != "1" or line[_MT].strip() != "451":
            continue
        try:
            za = int(round(endf_float(line[_FIELD[0]])))
            # Second record, fourth field: LIS0, which is 0 for a ground state.
            state = head[index + 1][_FIELD[3]].strip()
            if int(state or 0) != 0:
                return None
            match = _ZSYMAM_RE.match(head[index + _ZSYMAM][_FIELD[0]])
        except (ValueError, IndexError):
            return None
        if match is None:
            return None
        z, symbol, mass = int(match.group(1)), match.group(2).capitalize(), int(match.group(3))
        # ZSYMAM is free text, so check it against the header's own ZA.
        if z * 1000 + mass != za:
            return None
        return f"{symbol}{mass}"
    return None


def scan(endf_dir, wanted):
    """Find each wanted nuclide in a directory of ENDF files of any naming."""
    found = {}
    scanned = 0
    for path in sorted(endf_dir.rglob("*")):
        if not path.is_file():
            continue
        scanned += 1
        name = identify(path)
        if name is None or name not in wanted:
            continue
        if name in found:
            raise SystemExit(
                f"{name} is in {endf_dir} twice:\n  {found[name]}\n  {path}\n"
                "Point --endf-dir at one evaluation of each."
            )
        found[name] = path

    print(f"ENDF:  scanned {scanned} files in {endf_dir}")
    missing = [name for name in wanted if name not in found]
    if missing:
        raise SystemExit(
            f"no evaluation for {', '.join(missing)} in {endf_dir}"
            + (f" (found {', '.join(sorted(found))})" if found else "")
        )
    for name in wanted:
        print(f"       {name} <- {found[name].name}")
    return found


def tendl_name(nuclide):
    """"Fe56" -> "n-Fe056.tendl", TENDL's own file naming."""
    symbol = nuclide.rstrip("0123456789")
    mass = int(nuclide[len(symbol):])
    return f"n-{symbol}{mass:03d}.tendl"


def fetch_tendl(tarball, work_dir, nuclides, library):
    """The download fallback: unpack the wanted members of `library`'s archive.

    `work_dir` must already be scoped to `library`. TENDL names its evaluations
    identically in every release, so two releases unpacked into one directory
    would silently share files and the second library would never be read.

    The layout inside the archive differs between releases (TENDL-2025 is flat,
    TENDL-2017 nests under ``neutron_file/<El>/<El><A>/lib/endf/``), which costs
    nothing here because `_extract` matches on the basename.
    """
    wanted = {name: tendl_name(name) for name in nuclides}
    work_dir.mkdir(parents=True, exist_ok=True)

    have = {n: work_dir / f for n, f in wanted.items() if (work_dir / f).is_file()}
    if len(have) == len(wanted):
        print(f"ENDF:  already unpacked into {work_dir}")
        return have

    members = set(wanted.values())
    if tarball is not None and tarball.is_file():
        print(f"ENDF:  extracting {len(members)} files from {tarball}")
        with tarfile.open(tarball, "r:gz") as tar:
            _extract(tar, members, work_dir)
    else:
        url = data_source.tarball_url(library)
        print(f"ENDF:  streaming {url}")
        print(f"       (~3 GB read, only the {len(members)} wanted files are kept)")
        print("       Use --endf-dir to convert files you already have instead.")
        with urllib.request.urlopen(url) as response:
            with tarfile.open(fileobj=response, mode="r|gz") as tar:
                _extract(tar, members, work_dir)

    resolved = {n: work_dir / f for n, f in wanted.items()}
    for name, path in resolved.items():
        if not path.is_file():
            raise SystemExit(f"{name}: {path.name} was not in the archive")
    return resolved


def _is_whole_sublibrary(directory):
    """Whether a decay directory holds more than the metastable evaluations.

    Earlier versions of this script kept only the metastable files, which is
    all the isomer table reads. The reaction topology needs the whole
    sublibrary, so a directory left behind by one of those runs has to be
    filled in rather than trusted.
    """
    return any(not fnmatch.fnmatch(path.name, DECAY_ISOMER_GLOB)
               for path in directory.glob(DECAY_GLOB))


def fetch_decay(decay_dir, work_dir):
    """A directory holding the ENDF/B-8.1 decay sublibrary.

    Uses `decay_dir` if it already has one, otherwise pulls it. The isomer
    table reads only the 738 metastable evaluations, but the reaction topology
    needs every nuclide the network can reach, so all of them are kept.
    """
    if decay_dir is not None:
        decay_dir = pathlib.Path(decay_dir)
        if not _is_whole_sublibrary(decay_dir):
            raise SystemExit(f"--decay-dir {decay_dir} holds no {DECAY_GLOB} files, "
                             "or holds only the metastable ones")
        return decay_dir

    work_dir.mkdir(parents=True, exist_ok=True)
    if _is_whole_sublibrary(work_dir):
        return work_dir

    print(f"DECAY: fetching {DECAY_URL}")
    print("       (10.6 MB, 68 MB unpacked)")
    kept = 0
    with urllib.request.urlopen(DECAY_URL) as response:
        with tarfile.open(fileobj=response, mode="r|gz") as tar:
            for info in tar:
                name = pathlib.PurePosixPath(info.name).name
                if not info.isfile() or not fnmatch.fnmatch(name, DECAY_GLOB):
                    continue
                source = tar.extractfile(info)
                with open(work_dir / name, "wb") as out:
                    shutil.copyfileobj(source, out)
                kept += 1
    if not kept:
        raise SystemExit(f"no {DECAY_GLOB} files in {DECAY_URL}")
    print(f"       {kept} evaluations in {work_dir}")
    return work_dir


def build_branching(sources, decay_dir, out_root, library):
    """Write a `branching/` subsection from the same evaluations as the rates.

    Which fraction of an (n,2n) leaves the product in its metastable state
    rather than its ground state is energy dependent, and lives in MF=8/9/10 of
    the neutron evaluation. Without it the run falls back to whatever scalar
    branching the reference chain carries, which is wrong by up to 6x for foils
    whose decay heat comes from an isomer.

    Scoped to the foil's own isotopes. That covers the reactions that matter
    here, since one 5 minute irradiation barely burns the products, but it is
    not the whole chain.
    """
    print(f"BRANCH: isomeric branching for {len(sources)} parents")
    stats = yats.convert_branching(
        neutron_files=[str(path.resolve()) for _, path in sorted(sources.items())],
        decay_files=[str(p) for p in sorted(decay_dir.glob(DECAY_GLOB))],
        output_path=str(out_root),
        library=library,
    )
    interesting = {k: v for k, v in stats.items() if isinstance(v, int) and v}
    print(f"        {interesting} -> {out_root}")
    return out_root


def build_reactions(sources, decay_dir, out_root, library):
    """Write a `reactions/` subsection: the network topology, from your library.

    Which reaction on which parent gives which product, and the Q value of
    each, is something a neutron evaluation states directly, in the MT numbers
    it carries and their MF=3 headers. Taking it from the same files as the
    cross sections keeps the channel list and the rates consistent, and it is
    a much longer list: TENDL-2025 gives Fe56 25 reaction channels where the
    endf-b8.1 chain gives 5.

    Only the topology comes from the neutron files. Half-lives, decay modes and
    decay energies are read from `decay_dir` and stay on endf-b8.1, because no
    neutron-reaction library has them. Fission yields are not built at all;
    none of the FNS foils is fissionable.

    Scoped to the foil's own isotopes, for the same reason the branching is:
    one 5 minute irradiation at 1e10 n/cm2/s barely burns the products, so the
    parents that matter are the ones the foil started as.
    """
    print(f"REACT: reaction topology for {len(sources)} parents")
    yats.convert_transmutation(
        decay_files=[str(p) for p in sorted(decay_dir.glob(DECAY_GLOB))],
        fpy_files=[],
        neutron_files=[str(path.resolve()) for _, path in sorted(sources.items())],
        output_path=str(out_root),
        library=library,
        subsections=["reactions"],
    )
    print(f"        -> {out_root}")
    return out_root


def _extract(tar, members, dest):
    """Write just `members` out of an open tar, streaming-safe (no seeking)."""
    remaining = set(members)
    for info in tar:
        name = pathlib.PurePosixPath(info.name).name
        if name not in remaining:
            continue
        source = tar.extractfile(info)
        with open(dest / name, "wb") as out:
            shutil.copyfileobj(source, out)
        remaining.discard(name)
        print(f"       {name}")
        if not remaining:
            return


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--case", default="Fe",
                        help="FNS foil to convert for, e.g. Fe, W, SS316 (default: Fe)")
    parser.add_argument("--experiment", help="which experiment's deck to read the "
                                             "composition from (default: 2000exp_5min)")
    parser.add_argument("--fns-data", type=pathlib.Path, default=None,
                        help="benchmark JSON (default: ./fns_data.json)")
    parser.add_argument("--endf-dir", type=pathlib.Path,
                        help="directory of ENDF files, searched recursively; "
                             "also read from $ENDF_DIR")
    parser.add_argument("--library", default=data_source.DEFAULT_LIBRARY,
                        help="library name stamped into the Arrow output, and the "
                             "one downloaded without --endf-dir/--tarball. "
                             f"Downloadable: {', '.join(sorted(data_source.TENDL_TARBALLS))} "
                             f"(default: {data_source.DEFAULT_LIBRARY})")
    parser.add_argument("--output", type=pathlib.Path, default=None,
                        help="where the .arrow directories go "
                             "(default: data/<source>/neutron)")
    parser.add_argument("--chain", type=pathlib.Path, default=None,
                        help="where the branching and reaction subsections go "
                             "(default: data/<source>/chain)")
    parser.add_argument("--tarball", type=pathlib.Path,
                        help="local TENDL-n.tgz, used only without --endf-dir")
    parser.add_argument("--temperature", type=float, default=294.0, help="Kelvin (default: 294)")
    parser.add_argument("--njoy", default=shutil.which("njoy") or "njoy",
                        help="NJOY executable (default: njoy on PATH)")
    parser.add_argument("--decay-dir", type=pathlib.Path,
                        help="ENDF decay files for the isomer table; downloaded if absent")
    parser.add_argument("--no-branching", action="store_true",
                        help="skip the isomeric branching subsection")
    parser.add_argument("--no-reactions", action="store_true",
                        help="skip the reaction topology subsection")
    parser.add_argument("--force", action="store_true", help="reconvert files that already exist")
    parser.add_argument("--source", default=None,
                        help="folder name to file this run's data and results under "
                             "(default: the library name, or the --endf-dir directory name)")
    args = parser.parse_args()

    if shutil.which(args.njoy) is None and not pathlib.Path(args.njoy).is_file():
        raise SystemExit(
            f"NJOY not found ({args.njoy!r}). The njoy2016 wheel in requirements.txt "
            "puts one in the venv; pass --njoy to use a build of your own."
        )
    # njoy writes a listing called `output` into whatever directory it runs in,
    # so the version probe gets a temporary one to litter instead of this one.
    with tempfile.TemporaryDirectory() as scratch:
        banner = subprocess.run(
            [args.njoy], input="stop\n", capture_output=True, text=True, cwd=scratch
        ).stdout.splitlines()
    version = next((line.strip() for line in banner if line.strip()), "unknown version")
    print(f"NJOY:  {args.njoy} ({version})")

    case = fns_case.load(args.case, args.experiment, args.fns_data)
    print(f"CASE:  {case.describe()}")
    wanted = isotopes_for(case.elements)
    print(f"       {len(wanted)} natural isotopes to convert: {' '.join(wanted)}")

    endf_dir = args.endf_dir
    if endf_dir is None and os.environ.get("ENDF_DIR"):
        endf_dir = pathlib.Path(os.environ["ENDF_DIR"])

    # Everything this run writes goes under the name of the data it used, so a
    # second library cannot overwrite the first one's files (see data_source).
    source = data_source.slug(args.library, endf_dir, args.source)
    work_dir = HERE / "data" / source
    if args.output is None:
        args.output = work_dir / "neutron"
    if args.chain is None:
        args.chain = work_dir / "chain"
    print(f"SOURCE: {source} (data/{source}, results/{source})")

    if endf_dir is not None:
        if not endf_dir.is_dir():
            raise SystemExit(f"--endf-dir {endf_dir} is not a directory")
        sources = scan(endf_dir, wanted)
    else:
        sources = fetch_tendl(args.tarball, work_dir / "endf", wanted, args.library)

    out_dir = args.output
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"LIB:   stamping the output as {args.library!r}")

    for name, abundance in wanted.items():
        target = out_dir / f"{name}.arrow"
        if target.is_dir() and not args.force:
            print(f"{name}: already converted ({abundance:.4g}% of natural {name.rstrip('0123456789')})")
            continue
        print(f"{name}: NJOY at {args.temperature:g} K "
              f"({abundance:.4g}% of natural {name.rstrip('0123456789')}) ...", flush=True)
        start = time.perf_counter()
        yats.convert_neutron_xs(
            input_path=str(sources[name]),
            output_dir=str(out_dir),
            source_format="endf",
            temperatures=[args.temperature],
            library=args.library,
            njoy_exec=args.njoy,
        )
        size = sum(f.stat().st_size for f in target.rglob("*")) / 1e6
        print(f"{name}: {time.perf_counter() - start:.0f} s, {size:.1f} MB -> {target}")

    if args.no_branching and args.no_reactions:
        print("\nCHAIN: skipped, so isomer-dominated foils (Nb, Ag, W) will be wrong")
    else:
        decay_dir = fetch_decay(args.decay_dir, HERE / "data" / "decay")
        if args.no_branching:
            print("\nBRANCH: skipped, so isomer-dominated foils (Nb, Ag, W) will be wrong")
        else:
            build_branching(sources, decay_dir, args.chain, args.library)
        if args.no_reactions:
            print("REACT: skipped, so the network topology stays on " + DECAY_LIBRARY)
        else:
            build_reactions(sources, decay_dir, args.chain, args.library)

    print(f"\nArrow cross sections in {out_dir}")
    print(f"Next: python run_transmutation.py --case {args.case}"
          + (f" --experiment {args.experiment}" if args.experiment else ""))


if __name__ == "__main__":
    sys.exit(main())
