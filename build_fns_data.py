#!/usr/bin/env python3
"""Pack the FNS decay-heat benchmark into one JSON file.

The benchmark is distributed as three files per experiment across 73 foil
directories, 396 files in all, and most of that is the same four spectra
written out 132 times. This reads that tree once and writes `fns_data.json`,
which is what `fns_case.py` loads at runtime.

    build_fns_data.py --fns-dir ~/openmc_activator/fns

Run only when the benchmark data changes. The result is committed, so nobody
needs the original tree to use this example; the point of keeping this script
is that it says exactly how the JSON was derived, and lets anyone re-derive it.

Reading the FISPACT input decks happens here rather than at runtime. Only the
five keywords that describe the physics are taken from them (DENSITY, MASS, the
element lines, FLUX and TIME); the rest is solver settings and output requests,
which mean nothing outside FISPACT.
"""

import argparse
import collections
import json
import pathlib

# FISPACT time units. A bare TIME is in seconds, and a trailing ATOMS is an
# output request rather than a unit.
SECONDS = {
    "SECS": 1.0, "SECONDS": 1.0,
    "MINS": 60.0, "MINUTES": 60.0,
    "HOURS": 3600.0,
    "DAYS": 86400.0,
    "YEARS": 365.25 * 86400.0,
}

GROUPS = 709
GROUP_STRUCTURE = "CCFE-709"


def read_deck(path):
    """Density, mass, mass fractions and schedule from a FISPACT input deck."""
    density = mass_g = None
    composition = {}
    irradiation, cooling = [], []
    flux = 0.0
    expect_elements = 0

    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith(("*", "<<", "/")):
            continue
        tokens = line.split()
        keyword = tokens[0].upper()

        if expect_elements:
            # The lines right after MASS are `SYMBOL <weight percent>`.
            composition[tokens[0].capitalize()] = float(tokens[1]) / 100.0
            expect_elements -= 1
            continue

        if keyword == "DENSITY":
            density = float(tokens[1])
        elif keyword == "MASS":
            mass_g = float(tokens[1]) * 1000.0  # the deck is in kilograms
            expect_elements = int(tokens[2])
        elif keyword == "FLUX":
            flux = float(tokens[1])
        elif keyword == "TIME":
            duration = float(tokens[1])
            unit = tokens[2].upper() if len(tokens) > 2 else ""
            duration *= SECONDS.get(unit, 1.0)
            if flux > 0.0:
                irradiation.append([duration, flux])
            else:
                cooling.append(duration)

    if density is None or mass_g is None or not composition:
        raise SystemExit(f"{path}: no DENSITY/MASS/composition found")
    if not irradiation:
        raise SystemExit(f"{path}: no irradiation step found")
    total = sum(composition.values())
    if abs(total - 1.0) > 0.02:
        raise SystemExit(f"{path}: weight percentages sum to {total * 100:.2f}, not 100")
    return density, mass_g, composition, irradiation, cooling


def read_fluxes(path):
    """(709 group fluxes low energy first, irradiation position).

    The file is written highest group first, which is the order FISPACT wants
    and the reverse of the CCFE-709 boundaries, so the values are flipped here
    once rather than at every load. Its trailing comment names the position in
    the FNS assembly the foil sat at, which is what makes two foils share a
    spectrum.
    """
    values, position = [], None
    for line in path.read_text().splitlines():
        if "Position" in line:
            after = line.split("Position", 1)[1].strip().lstrip("#").split()
            position = after[0] if after else None
            break
        if "Normalization" in line:
            break
        values.extend(float(token) for token in line.split())
    if len(values) < GROUPS:
        raise SystemExit(f"{path}: expected {GROUPS} group fluxes, found {len(values)}")
    if position is None:
        raise SystemExit(f"{path}: no irradiation position in the trailing comment")
    return values[:GROUPS][::-1], position


def read_measurement(path):
    """(times, specific decay heat [uW/g], uncertainty) from a .exp file.

    Times are cumulative from the end of the irradiation, in a unit the file
    does not state. Some files are padded to a fixed length with all-zero rows,
    which are not measurements.
    """
    times, heat, uncertainty = [], [], []
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) != 3:
            continue
        time, value, sigma = (float(p) for p in parts)
        if time == 0.0 and value == 0.0:
            continue
        times.append(time)
        heat.append(value)
        uncertainty.append(sigma)
    if not times:
        raise SystemExit(f"{path}: no measured points")
    return times, heat, uncertainty


def time_unit(deck_cooling, last_time, path):
    """Which unit the .exp times are in, judged against the deck's schedule.

    The measurement file gives bare numbers. The deck's cooling steps cover
    roughly the same span in units it does state, so whichever unit puts the
    last measured point closest to the end of that span is the one meant. The
    candidates are 60x apart, so a deck that disagrees with its measurement in
    detail still picks the right one. This is the only thing the deck's own
    cooling schedule is needed for.
    """
    span = sum(deck_cooling)
    if span <= 0.0 or last_time <= 0.0:
        raise SystemExit(f"{path}: cannot tell what unit the measured times are in")
    candidates = {"minutes": 60.0, "hours": 3600.0, "days": 86400.0}
    name = min(candidates, key=lambda u: abs(last_time * candidates[u] / span - 1.0))
    ratio = last_time * candidates[name] / span
    if not 0.5 < ratio < 2.0:
        raise SystemExit(
            f"{path}: measured times end at {last_time:g} {name}, but the deck's "
            f"cooling covers {span:g} s. These do not describe the same experiment."
        )
    return name


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fns-dir", type=pathlib.Path,
                        default=pathlib.Path.home() / "openmc_activator" / "fns",
                        help="the benchmark tree, one directory per foil")
    parser.add_argument("-o", "--output", type=pathlib.Path,
                        default=pathlib.Path(__file__).resolve().parent / "fns_data.json")
    args = parser.parse_args()

    if not args.fns_dir.is_dir():
        raise SystemExit(f"{args.fns_dir} is not a directory")

    spectra, cases = {}, {}
    experiments = collections.Counter()
    for directory in sorted(d for d in args.fns_dir.iterdir() if d.is_dir()):
        for deck in sorted(directory.glob("*.i")):
            # TENDL-2017_1996exp_5min.i -> 1996exp_5min
            name = deck.stem.split("_", 1)[-1]
            fluxes = directory / f"{name}_fluxes"
            measurement = directory / f"{name}.exp"
            if not (fluxes.is_file() and measurement.is_file()):
                continue

            density, mass_g, composition, irradiation, cooling = read_deck(deck)
            groups, position = read_fluxes(fluxes)
            times, measured, uncertainty = read_measurement(measurement)

            # The four distinct spectra are stored once and referred to by the
            # position they were measured at, rather than repeated 132 times.
            if position in spectra:
                if spectra[position] != groups:
                    raise SystemExit(
                        f"{fluxes}: position {position} already has a different spectrum"
                    )
            else:
                spectra[position] = groups

            cases.setdefault(directory.name, {})[name] = {
                "density_g_per_cm3": density,
                "mass_g": mass_g,
                "composition_mass_fraction": composition,
                "irradiation_s_and_flux": irradiation,
                "position": position,
                "time_unit": time_unit(cooling, times[-1], measurement),
                "times": times,
                "measured_uW_per_g": measured,
                "uncertainty_uW_per_g": uncertainty,
            }
            experiments[name] += 1

    document = {
        "description": "JAEA FNS decay-heat benchmark, as distributed through the "
                       "IAEA CoNDERC fusion benchmark collection.",
        "source": "https://github.com/jbae11/openmc_activator (fns/)",
        "built_by": "build_fns_data.py",
        "group_structure": GROUP_STRUCTURE,
        "spectrum_note": f"{GROUPS} absolute group fluxes [n/cm2/s], lowest energy "
                         f"first, matching the {GROUP_STRUCTURE} boundaries. Their sum "
                         "is the total flux the foil saw.",
        "spectra_by_position": spectra,
        "cases": cases,
    }
    args.output.write_text(json.dumps(document, indent=1, sort_keys=True))

    size = args.output.stat().st_size / 1024
    print(f"{args.output}: {size:.0f} KB")
    print(f"  {len(cases)} foils, {sum(experiments.values())} experiments, "
          f"{len(spectra)} distinct spectra")
    for name, count in sorted(experiments.items()):
        print(f"  {name:16s} {count:3d} foils")


if __name__ == "__main__":
    main()
