#!/usr/bin/env python3
"""Step 2: irradiate one FNS foil with yani and plot against the measurement.

The FNS decay-heat experiments (JAEA Fusion Neutronics Source, IAEA CoNDERC)
hold a 1 g foil in a 14 MeV field, pull it out, and measure how much heat it
gives off as the activation products decay.

    run_transmutation.py                                 # iron, the default
    run_transmutation.py --case W
    run_transmutation.py --case SS316
    run_transmutation.py --case Fe --experiment 1996exp_7hour

Everything about the case (composition, density, flux, spectrum, schedule and
the measurement) comes from `fns_data.json`. Cross sections come from the Arrow
directory convert_to_arrow.py produced, which has to hold the same foil's
isotopes, as do the reaction topology and isomeric branching. Decay data comes
from endf-b8.1, which yani downloads on first use.

There is no transport here: the measured spectrum is the input, and yani
collapses it against the cross sections to get one-group reaction rates, then
solves the Bateman system over the schedule.
"""

import argparse
import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import yani  # noqa: E402

import data_source  # noqa: E402
import fns_case  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent

# The temperature the cross sections were broadened to by convert_to_arrow.py.
TEMPERATURE = 294.0

# Half-lives, decay modes and decay energies, and the fallback for anything
# convert_to_arrow.py did not write. Held on endf-b8.1 because TENDL publishes
# no decay sublibrary.
DECAY_LIBRARY = "endf-b8.1"


def available_cross_section_sources(data_root):
    """Source names that already have a converted neutron Arrow directory."""
    if not data_root.is_dir():
        return []
    return sorted(path.parent.name for path in data_root.glob("*/neutron") if path.is_dir())


def inferred_chain_path(cross_sections):
    """The matching chain directory for a standard neutron path, or None.

    convert_to_arrow.py writes cross sections to data/<source>/neutron and the
    chain subsections to data/<source>/chain. If run_transmutation.py is pointed
    at that neutron directory, prefer the sibling chain by default so both come
    from the same conversion run.
    """
    if cross_sections is None:
        return None
    if pathlib.Path(cross_sections).name == "neutron":
        return pathlib.Path(cross_sections).parent / "chain"
    return None


def subsections_of(chain):
    """Which chain subsections convert_to_arrow.py wrote, from the manifest."""
    if chain is None or not (chain / "manifest.json").is_file():
        return set()
    return set(json.loads((chain / "manifest.json").read_text()).get("subsections", {}))


def scope_of(chain):
    """Which parents the chain covers, from the sidecar convert_to_arrow.py wrote.

    None where the chain predates the sidecar, which is not the same as a chain
    covering nothing: those fall through to the zero-heat check after the solve.
    """
    if chain is None or not (chain / "scope.json").is_file():
        return None
    return set(json.loads((chain / "scope.json").read_text()).get("parents", []))


def run(case, cross_sections, chain):
    """Specific decay heat [uW/g] after each cooling step, and its breakdown."""
    yani.cross_section_data = str(cross_sections)
    yani.transmutation_decay_data = DECAY_LIBRARY
    yani.transmutation_fission_yields = DECAY_LIBRARY

    # Everything a neutron evaluation states is taken from the one the cross
    # sections came from, so that the whole reaction side of the network moves
    # together when the library changes. That is the topology (which reaction
    # on which parent gives which product) and the isomeric branching (which
    # state it lands in, which is energy dependent and decides the answer
    # outright for foils whose heat comes from an isomer).
    have = subsections_of(chain)
    for subsection, setting, warning in [
        ("branching", "transmutation_branch_ratios",
         "isomer-dominated foils will be off; rerun without --no-branching"),
        ("reactions", "transmutation_reactions",
         "fewer production channels; rerun without --no-reactions"),
    ]:
        if subsection in have:
            setattr(yani, setting, str(chain))
            print(f"{subsection}: {chain}")
        else:
            setattr(yani, setting, DECAY_LIBRARY)
            print(f"{subsection}: {DECAY_LIBRARY} ({warning})")

    # The histogram is a shape; the pulse rate carries the magnitude.
    spectrum = yani.NeutronSource(
        energy=yani.sources.Histogram("CCFE-709", case.spectrum.tolist())
    )
    schedule = yani.PulseSchedule(
        [yani.Pulse(duration=(seconds, "s"), rate=flux, source=spectrum)
         for seconds, flux in case.irradiation]
        + [yani.Cooldown(duration=(seconds, "s")) for seconds in case.cooling]
    )

    material = yani.Material(
        composition=case.composition,
        fraction_type="mass",
        density=case.density,
        volume=case.mass_g / case.density,
        temperature=TEMPERATURE,
        name=f"{case.name} foil",
    )
    nuclides = sorted(material.get_nuclide_names())
    print(f"material: {len(nuclides)} nuclides, {' '.join(nuclides)}")

    # The chain is rebuilt per foil and scoped to that foil's isotopes, so one
    # left over from a different foil carries no reactions for these nuclides
    # and solves to an inventory of nothing. That surfaced only as zero decay
    # heat at the end of a full run, which names neither the foil the chain was
    # built for nor the nuclides it is missing. Refuse it up front instead.
    covered = scope_of(chain)
    missing = sorted(set(nuclides) - covered) if covered is not None else []
    if missing:
        raise SystemExit(
            f"{case.name}: the chain at {chain} was built for "
            f"{' '.join(sorted(covered)) or 'nothing'}, so it has no reactions for "
            f"{' '.join(missing)}.\n"
            f"Rebuild it for this foil and try again:\n"
            f"  python convert_to_arrow.py --case {case.name} "
            f"--output {cross_sections} --chain {chain}"
        )

    # `transmute` hands back a TransmutationResults keyed by material id, not a
    # plain list. Its index 0 is the initial composition, so `step_materials` is
    # the per-step view the slice below counts from. The foil is built without an
    # explicit id, which files it under 0.
    results = material.transmute(schedule)
    states = results.step_materials(material.id or 0)

    # One state per schedule step and no initial entry, so the irradiation
    # pulses come first and the measurement starts after them.
    heat, breakdown = [], []
    for state in states[len(case.irradiation):]:
        contributions = state.decay_heat(by_nuclide=True)
        heat.append(sum(contributions.values()) / case.mass_g * 1e6)
        breakdown.append(contributions)
    heat = np.array(heat)

    # The rate of every production edge the solve drove, which yani-core 0.9.0
    # hands back and earlier versions computed and threw away
    # (fusion-neutronics/core#505). Summed over the irradiation pulses and
    # weighted by their duration, so an edge carries the production it drove per
    # atom of its parent over the whole irradiation, which is what weights a
    # route. Cooldown steps drive no reactions and contribute nothing.
    #
    # Ratios within one (parent, kind) are the flux-weighted isomeric branching:
    # the split arrives as two edges of one channel, so the sum is the channel
    # rate and the ratio is f_m. That is the number the chain file cannot give,
    # because its own branching for the dominant tungsten channel is a
    # placeholder the energy-dependent overlay replaces at solve time.
    edges = {}
    for step, (seconds, _flux) in enumerate(case.irradiation):
        for parent, kinds in (results.get_reaction_rates(material.id or 0, step)
                              or {}).items():
            for kind, targets in kinds.items():
                for target, rate in targets:
                    key = (parent, kind, target)
                    edges[key] = edges.get(key, 0.0) + rate * seconds
    edge_rates = [[parent, kind, target, weight]
                  for (parent, kind, target), weight in sorted(
                      edges.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2] or ""))]

    # The parents a route starts from, needed to turn a per-parent-atom edge
    # weight into a share of what actually reached a product. Index 0 of the
    # results is the composition before any step; `states` deliberately skips
    # it, so this asks the results rather than reusing `states[0]`, which is the
    # material after the first pulse and already burnt.
    initial = results.get_material(material.id or 0, 0)
    initial_atoms = dict(initial.nuclides) if initial is not None else {}

    # An irradiated foil always activates, so heat that is identically zero
    # means the network produced nothing at all. It used to be reported as
    # C/E 0.000, which reads like a physics result rather than the wiring
    # mistake it invariably is, and a sweep would log it foil after foil
    # without ever failing. It is an error here instead.
    #
    # The chain is the usual suspect: it is rebuilt per foil and scoped to that
    # foil's isotopes, so one built for a different foil leaves these nuclides
    # with no reactions. A chain from an older converter, or cross sections that
    # do not cover the material, will do it too.
    if not heat.any():
        raise SystemExit(
            f"{case.name}: every cooling step came out at zero decay heat, so the "
            f"transmutation produced no activation products.\n"
            f"  chain:          {chain}\n"
            f"  cross sections: {cross_sections}\n"
            f"Rebuild both for this foil and try again:\n"
            f"  python convert_to_arrow.py --case {case.name} "
            f"--output {cross_sections} --chain {chain}"
        )
    return heat, breakdown, edge_rates, initial_atoms


def summarise(case, calculated):
    """C/E per point, plus the one-line numbers that go in the plot title."""
    ratio = np.divide(calculated, case.measured,
                      out=np.full_like(calculated, np.nan), where=case.measured > 0)
    sigma = case.uncertainty / np.where(case.measured > 0, case.measured, np.nan)
    return ratio, {
        "median_ratio": float(np.nanmedian(ratio)),
        "mean_deviation_percent": float(np.nanmean(np.abs(ratio - 1.0)) * 100.0),
        "median_measurement_sigma_percent": float(np.nanmedian(sigma) * 100.0),
    }


def plot(case, calculated, breakdown, ratio, metrics, out_dir):
    times, measured, uncertainty = case.times, case.measured, case.uncertainty
    figure, (top, bottom) = plt.subplots(
        2, 1, figsize=(8, 8), sharex=True, gridspec_kw={"height_ratios": [3, 1]}
    )

    handles = [
        top.errorbar(times, measured, yerr=uncertainty, fmt="o", color="black",
                     markersize=5, capsize=3, label="FNS measurement", zorder=3),
        top.plot(times, calculated, "-", color="#d62728", linewidth=2,
                 label="yani", zorder=2)[0],
    ]

    # Name the products that carry the heat, so the curve is readable as physics
    # rather than as one number.
    totals = {}
    for contributions in breakdown:
        for nuclide, watts in contributions.items():
            totals[nuclide] = totals.get(nuclide, 0.0) + watts
    # Only the leading few are drawn, to keep the legend readable. The rest are
    # not lost: every nuclide the network produced is written to the JSON under
    # by_nuclide_uW_per_g, one entry per cooling step, so anything dropped here
    # can still be plotted or ranked from the saved result.
    #
    # The ranking is by heat summed over the cooling steps, which favours
    # whatever dominates the early high-heat points. A product that peaks late
    # can carry a visible share of the curve and still fall outside the cut, so
    # this legend is not the list of everything that matters.
    leaders = sorted(totals, key=totals.get, reverse=True)[:3]
    for nuclide in leaders:
        series = [c.get(nuclide, 0.0) / case.mass_g * 1e6 for c in breakdown]
        handles.append(top.plot(times, series, "--", linewidth=1, alpha=0.7, label=nuclide)[0])

    top.set_ylabel("specific decay heat [uW/g]")
    top.set_yscale("log")
    # Short-lived contributors fall away by orders of magnitude; without a floor
    # they stretch the axis until the total and the measurement are one line.
    positive = calculated[calculated > 0]
    if positive.size:
        top.set_ylim(bottom=positive.min() / 300.0, top=positive.max() * 3.0)
    # The agreement goes on the figure, not just on stdout, so a directory of
    # these can be read side by side. The measurement's own scatter sits next to
    # the deviation because it says whether that deviation means anything.
    irradiation = sum(seconds for seconds, _ in case.irradiation)
    top.set_title(
        f"FNS {case.name}, {case.experiment}: "
        f"{irradiation / 60:g} minute irradiation at 14 MeV\n"
        f"median C/E {metrics['median_ratio']:.3f}, "
        f"mean deviation {metrics['mean_deviation_percent']:.1f}%, "
        f"measurement sigma {metrics['median_measurement_sigma_percent']:.1f}%",
        fontsize=11,
    )
    top.legend(handles=handles, loc="lower left")
    top.grid(alpha=0.3)

    bottom.axhline(1.0, color="black", linewidth=1)
    bottom.fill_between(times, 1 - uncertainty / measured, 1 + uncertainty / measured,
                        color="black", alpha=0.15, label="measurement uncertainty")
    bottom.plot(times, ratio, "o-", color="#d62728", markersize=4)
    bottom.set_xlabel(f"time after shutdown [{case.time_unit}]")
    bottom.set_ylabel("yani / measured")
    bottom.set_xscale("log")
    bottom.legend()
    bottom.grid(alpha=0.3)

    figure.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"fns_{case.name}_{case.experiment}.png"
    figure.savefig(path, dpi=150)
    plt.close(figure)
    print(f"\nplot: {path}")
    return ratio


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--case", default="Fe",
                        help="FNS foil, e.g. Fe, W, SS316 (default: Fe)")
    parser.add_argument("--experiment", help="which experiment (default: 2000exp_5min)")
    parser.add_argument("--fns-data", type=pathlib.Path, default=None,
                        help="benchmark JSON (default: ./fns_data.json)")
    parser.add_argument("--library", default=data_source.DEFAULT_LIBRARY,
                        help="which data build to read and where to write results "
                             f"(default: {data_source.DEFAULT_LIBRARY})")
    parser.add_argument("--endf-dir", type=pathlib.Path,
                        help="names the build made from your own evaluations, so the "
                             "defaults below point at it")
    parser.add_argument("--cross-sections", type=pathlib.Path, default=None,
                        help="Arrow directory from convert_to_arrow.py "
                             "(default: data/<source>/neutron)")
    parser.add_argument("--chain", type=pathlib.Path, default=None,
                        help="branching subsection from convert_to_arrow.py "
                            "(default: data/<source>/chain, or the sibling of "
                            "--cross-sections when it ends with /neutron)")
    parser.add_argument("--output", type=pathlib.Path, default=None,
                        help="where the results go (default: results/<source>)")
    parser.add_argument("--list", action="store_true",
                        help="list the available foils and experiments, then exit")
    parser.add_argument("--source", default=None,
                        help="folder name to file this run's data and results under "
                             "(default: the library name, or the --endf-dir directory name)")
    args = parser.parse_args()

    # Results and inputs live under the name of the nuclear data they came from,
    # so two libraries do not overwrite each other (see data_source).
    source = data_source.slug(args.library, args.endf_dir, args.source)
    if args.cross_sections is None:
        args.cross_sections = HERE / "data" / source / "neutron"
    if args.chain is None:
        args.chain = inferred_chain_path(args.cross_sections)
        if args.chain is None:
            args.chain = HERE / "data" / source / "chain"
    if args.output is None:
        args.output = HERE / "results" / source

    if args.list:
        for name in fns_case.cases(args.fns_data):
            print(f"{name:8s} {' '.join(fns_case.experiments(name, args.fns_data))}")
        return

    if not args.cross_sections.is_dir():
        data_root = HERE / "data"
        available = available_cross_section_sources(data_root)
        message = [
            f"no cross sections at {args.cross_sections}.",
            "Run convert_to_arrow.py first, or point this run at an existing source.",
        ]
        if available:
            message.append(
                "Known converted sources under data/: " + ", ".join(available)
            )
            if source not in available:
                pick = available[0] if len(available) == 1 else None
                if pick is not None:
                    message.extend([
                        "Your source defaults to the library name unless --source is set.",
                        "This often differs after --endf-dir, where the directory basename "
                        "becomes the source name.",
                        "Try one of:",
                        f"  python run_transmutation.py --case {args.case} --source {pick}",
                        f"  python run_transmutation.py --case {args.case} "
                        f"--cross-sections data/{pick}/neutron",
                    ])
                else:
                    message.extend([
                        "Your source defaults to the library name unless --source is set.",
                        "This often differs after --endf-dir, where the directory basename "
                        "becomes the source name.",
                        "Use one of the listed source names with --source, for example:",
                        f"  python run_transmutation.py --case {args.case} --source <one-of-above>",
                    ])
        raise SystemExit("\n".join(message))

    case = fns_case.load(args.case, args.experiment, args.fns_data)
    print(f"case: {case.describe()}")
    print(f"spectrum: {fns_case.GROUPS} CCFE groups, {case.spectrum.sum():.4g} n/cm2/s summed")

    calculated, breakdown, edge_rates, initial_atoms = run(
        case, args.cross_sections, args.chain)
    ratio, metrics = summarise(case, calculated)
    plot(case, calculated, breakdown, ratio, metrics, args.output)

    print(f"\n{case.time_unit:>9}  {'measured':>10}  {'yani':>10}  {'C/E':>6}")
    for time, exp, calc, r in zip(case.times, case.measured, calculated, ratio):
        print(f"{time:9.2f}  {exp:10.3e}  {calc:10.3e}  {r:6.3f}")
    print(f"\nmedian C/E {metrics['median_ratio']:.3f}, "
          f"mean deviation {metrics['mean_deviation_percent']:.1f}%, "
          f"measurement sigma {metrics['median_measurement_sigma_percent']:.1f}%")

    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / f"fns_{case.name}_{case.experiment}.json").write_text(json.dumps({
        "case": case.name,
        "experiment": case.experiment,
        "time_unit": case.time_unit,
        "times": case.times.tolist(),
        "measured_uW_per_g": case.measured.tolist(),
        "measured_uncertainty": case.uncertainty.tolist(),
        "yani_uW_per_g": calculated.tolist(),
        "ratio": ratio.tolist(),
        **metrics,
        "by_nuclide_uW_per_g": [
            {n: w / case.mass_g * 1e6 for n, w in step.items()} for step in breakdown
        ],
        # [parent, reaction kind, target (null for fission), production per atom
        # of the parent over the whole irradiation]. make_report.py weights the
        # production routes with these; see `run` for what the number is.
        "edge_rates": edge_rates,
        "initial_atoms_per_barn_cm": initial_atoms,
    }, indent=2))


if __name__ == "__main__":
    main()
