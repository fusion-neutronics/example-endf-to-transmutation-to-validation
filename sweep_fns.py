#!/usr/bin/env python3
"""Step 3: every FNS foil, one library, one table.

Runs convert_to_arrow.py and run_transmutation.py over all 73 foils and
collects what came out into a ranked table. One foil says whether a handful of
cross sections are right; 73 of them say whether a library is.

    sweep_fns.py --endf-dir /path/to/tendl-2025
    sweep_fns.py --endf-dir /path/to/endfb-8.1 --library endf-b8.1

A local copy of the library is required, either `--endf-dir`, `$ENDF_DIR` or
`--tarball`. convert_to_arrow.py on its own will stream TENDL-n.tgz for the few
isotopes one foil needs, which is a reasonable trade once and a terrible one
73 times.

Both scripts are called as they would be by hand, so a foil that fails is
reported and the sweep carries on. Conversion runs each time, so a second run
rebuilds every foil from the same command line.

The table is sorted by mean deviation, worst last, because the point of running
all of them is to find where the library falls over. Read it alongside
`median_measurement_sigma_percent`: a foil measured to 30% cannot tell you much,
and there are a few of those.
"""

import argparse
import json
import os
import pathlib
import subprocess
import sys
import time

import data_source
import fns_case

HERE = pathlib.Path(__file__).resolve().parent


def run_case(case, experiment, results, converter_args, runner_args):
    """Convert and run one foil. Returns its summary, or None if it failed."""
    for script, extra in [("convert_to_arrow.py", converter_args),
                          ("run_transmutation.py", runner_args)]:
        command = [sys.executable, str(HERE / script), "--case", case,
                   "--experiment", experiment, *extra]
        finished = subprocess.run(command, capture_output=True, text=True)
        if finished.returncode != 0:
            tail = (finished.stderr or finished.stdout).strip().splitlines()
            print(f"  FAILED in {script}: {tail[-1] if tail else 'no output'}")
            return None

    written = results / f"fns_{case}_{experiment}.json"
    if not written.is_file():
        print(f"  FAILED: {written} was not written")
        return None
    return json.loads(written.read_text())


def data_sigma(summary):
    """The median nuclear-data sigma on the calculated heat, as a percent, or None.

    Beside the measurement's own sigma this is what says whether a foil's
    deviation is a result. A 10% deviation on a foil whose cross sections are
    known to 3% is a disagreement; the same 10% on one known to 25% is not, and
    the sweep table cannot rank foils honestly without both numbers.

    ``None`` for a foil run with ``--no-uncertainty``, or against cross sections
    converted without the covariance.
    """
    sigma = summary.get("yani_uncertainty_uW_per_g")
    if sigma is None:
        return None
    calculated = summary["yani_uW_per_g"]
    relative = [100.0 * s / c for s, c in zip(sigma, calculated) if c > 0]
    if not relative:
        return None
    relative.sort()
    middle = len(relative) // 2
    return (relative[middle] if len(relative) % 2
            else 0.5 * (relative[middle - 1] + relative[middle]))


def dominant(summary):
    """The nuclide carrying most of the heat at the first measured point."""
    contributions = summary["by_nuclide_uW_per_g"][0]
    total = sum(contributions.values())
    if total <= 0.0:
        return "-", 0.0
    nuclide, heat = max(contributions.items(), key=lambda item: item[1])
    return nuclide, heat / total * 100.0


def table(rows, library, experiment):
    """The sweep as markdown, worst agreement last."""
    lines = [
        "| foil | median C/E | mean deviation | measurement sigma | data sigma "
        "| dominant product |",
        "|---|---|---|---|---|---|",
    ]
    for row in sorted(rows, key=lambda r: r["mean_deviation_percent"]):
        sigma = row["median_data_sigma_percent"]
        lines.append(
            f"| {row['case']} | {row['median_ratio']:.3f} | "
            f"{row['mean_deviation_percent']:.1f}% | "
            f"{row['median_measurement_sigma_percent']:.1f}% | "
            f"{'-' if sigma is None else f'{sigma:.1f}%'} | "
            f"{row['dominant']} {row['dominant_percent']:.0f}% |"
        )
    return (f"# FNS decay heat, {library} on {experiment}\n\n"
            f"{len(rows)} foils.\n\n" + "\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--experiment", default="2000exp_5min",
                        help="which campaign to run (default: 2000exp_5min, the only "
                             "one covering all 73 foils)")
    parser.add_argument("--cases", nargs="+", help="only these foils (default: all of them)")
    parser.add_argument("--endf-dir", type=pathlib.Path,
                        help="passed to convert_to_arrow.py; also read from $ENDF_DIR")
    parser.add_argument("--tarball", type=pathlib.Path,
                        help="local TENDL-n.tgz, used instead of --endf-dir")
    parser.add_argument("--library", default=data_source.DEFAULT_LIBRARY,
                        help="passed to convert_to_arrow.py, and the name this "
                             "sweep's data and results are filed under. "
                             f"Downloadable: {', '.join(sorted(data_source.TENDL_TARBALLS))} "
                             f"(default: {data_source.DEFAULT_LIBRARY})")
    parser.add_argument("--cross-sections", type=pathlib.Path, default=None,
                        help="Arrow directory (default: data/<source>/neutron)")
    parser.add_argument("--output", type=pathlib.Path, default=None,
                        help="where the per-foil results and the table go "
                             "(default: results/<source>/sweep)")
    parser.add_argument("--source", default=None,
                        help="folder name to file this run's data and results under "
                             "(default: the library name, or the --endf-dir directory name)")
    args = parser.parse_args()

    # One folder per data source, so a TENDL-2017 sweep and a TENDL-2025 sweep
    # sit side by side and can be compared (see data_source).
    source = data_source.slug(args.library, args.endf_dir, args.source)
    if args.cross_sections is None:
        args.cross_sections = HERE / "data" / source / "neutron"
    if args.output is None:
        args.output = HERE / "results" / source / "sweep"
    print(f"SOURCE: {source} -> {args.output}")

    wanted = args.cases or fns_case.cases()
    missing = [c for c in wanted if args.experiment not in fns_case.experiments(c)]
    if missing:
        raise SystemExit(f"{args.experiment} does not cover {', '.join(missing)}")

    if args.endf_dir is None and args.tarball is None and not os.environ.get("ENDF_DIR"):
        raise SystemExit(
            "the sweep needs a local copy of the library: --endf-dir, $ENDF_DIR or "
            "--tarball.\nWithout one, convert_to_arrow.py streams TENDL's 3.5 GB "
            f"archive once for each of the {len(wanted)} foils in this sweep."
        )

    # The chain is rebuilt per foil, scoped to that foil's own isotopes, so the
    # writer and the reader have to be handed the SAME path explicitly. Leaving
    # either of them to work it out from its own arguments is how a foil ends up
    # read against the PREVIOUS foil's chain: no reaction in it has this foil's
    # nuclides as a parent, so nothing is produced and the decay heat comes out
    # a silent zero rather than an error.
    chain_dir = HERE / "data" / source / "chain"
    converter_args = ["--library", args.library,
                      "--source", source,
                      "--output", str(args.cross_sections),
                      "--chain", str(chain_dir)]
    if args.endf_dir is not None:
        converter_args += ["--endf-dir", str(args.endf_dir)]
    if args.tarball is not None:
        converter_args += ["--tarball", str(args.tarball)]
    runner_args = ["--cross-sections", str(args.cross_sections),
                   "--chain", str(chain_dir),
                   "--output", str(args.output)]
    args.output.mkdir(parents=True, exist_ok=True)

    rows = []
    failed = []
    start = time.perf_counter()
    for index, case in enumerate(wanted, 1):
        elapsed = time.perf_counter() - start
        print(f"[{index}/{len(wanted)}] {case} ({elapsed / 60:.0f} min in)", flush=True)
        summary = run_case(case, args.experiment, args.output, converter_args, runner_args)
        if summary is None:
            failed.append(case)
            continue
        nuclide, percent = dominant(summary)
        sigma = data_sigma(summary)
        rows.append({
            "case": case,
            "median_ratio": summary["median_ratio"],
            "mean_deviation_percent": summary["mean_deviation_percent"],
            "median_measurement_sigma_percent": summary["median_measurement_sigma_percent"],
            "median_data_sigma_percent": sigma,
            "dominant": nuclide,
            "dominant_percent": percent,
        })
        print(f"  C/E {summary['median_ratio']:.3f}, "
              f"{summary['mean_deviation_percent']:.1f}% deviation"
              f"{'' if sigma is None else f' against {sigma:.1f}% data sigma'}, "
              f"{nuclide} {percent:.0f}%",
              flush=True)

    (args.output / "sweep.json").write_text(json.dumps(
        {"library": args.library, "experiment": args.experiment,
         "failed": failed, "foils": rows}, indent=1, sort_keys=True))
    (args.output / "sweep.md").write_text(table(rows, args.library, args.experiment))

    print(f"\n{len(rows)} foils in {(time.perf_counter() - start) / 60:.0f} min "
          f"-> {args.output}/sweep.md")
    if failed:
        print(f"{len(failed)} failed: {' '.join(failed)}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
