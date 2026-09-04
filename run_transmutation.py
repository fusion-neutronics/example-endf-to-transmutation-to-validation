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

The decay heat comes back with an uncertainty on it. yani resamples the
activation cross sections from the MF=33 covariance that step 1 wrote beside
each nuclide, folds each draw against this foil's own spectrum and re-solves the
schedule, so the spread over the ensemble is the nuclear-data uncertainty on the
answer. ``--no-uncertainty`` skips it and costs nothing to omit; the mean
inventories are identical either way.
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

# Peak share of the decay heat a product needs before its production routes are
# worked out and filed. Set under make_report.py's own SHARE_FLOOR so the report
# still decides what to print, and above zero because a chain reaches hundreds of
# nuclides and almost none of them ever carry a row.
ROUTE_PRODUCT_FLOOR = 0.001

# How far a route may run: neutron reactions, and then decay steps after them.
#
# One reaction is what the published pathway pages carry, and what these
# irradiations can support. A second reaction step costs another factor of the
# fluence, and at 1e10 n/cm2/s for five minutes that is a fluence of 3e12, so a
# two-step route arrives at parts in a billion of a one-step one.
#
# Three decay steps covers what these reach: the longest route the published
# pages carry is W182(n,p)Ta182n(IT)Ta182m(IT)Ta182, which is three.
ROUTE_REACTION_STEPS = 1
ROUTE_DECAY_STEPS = 3

# Half-lives, decay modes and decay energies, and the fallback for anything
# convert_to_arrow.py did not write. Held on endf-b8.1 because TENDL publishes
# no decay sublibrary.
DECAY_LIBRARY = "endf-b8.1"


def available_cross_section_sources(data_root):
    """Source names that already have a converted neutron Arrow directory."""
    if not data_root.is_dir():
        return []
    return sorted(path.parent.name for path in data_root.glob("*/neutron") if path.is_dir())


def inferred_chain_path(cross_sections, case):
    """The matching chain directory for a standard neutron path, or None.

    convert_to_arrow.py writes cross sections to data/<source>/neutron and the
    chain subsections to data/<source>/chain-<case>. If run_transmutation.py is
    pointed at that neutron directory, prefer the sibling chain by default so
    both come from the same conversion run.

    The foil is part of the name because the chain is scoped to one foil's
    isotopes while the cross sections are not: one library's neutron directory
    holds every nuclide ever converted into it, and one chain beside it can
    only ever be the last foil's.
    """
    if cross_sections is None:
        return None
    if pathlib.Path(cross_sections).name == "neutron":
        return pathlib.Path(cross_sections).parent / f"chain-{case}"
    return None


def subsections_of(chain):
    """Which chain subsections convert_to_arrow.py wrote, from the manifest."""
    if chain is None or not (chain / "manifest.json").is_file():
        return set()
    return set(json.loads((chain / "manifest.json").read_text()).get("subsections", {}))


def library_of(root, nuclides=None):
    """The library name the converted data on disk stamps itself with.

    Read out of the data rather than taken from the directory name. `--source`
    names a folder, and nothing stops a folder called one thing from holding
    evaluations converted from another; a report that says which library a
    number came from should be quoting the data and not the path.

    A chain says so in its manifest. A neutron directory says so once per
    nuclide, and only the foil's own nuclides are asked, because that directory
    accumulates every nuclide ever converted into it and the ones this run did
    not read have no bearing on this run.
    """
    root = pathlib.Path(root)
    manifest = root / "manifest.json"
    if manifest.is_file():
        return json.loads(manifest.read_text()).get("library") or None
    stamped = set()
    for nuclide in nuclides or []:
        version = root / f"{nuclide}.arrow" / "version.json"
        if version.is_file():
            stamped.add(json.loads(version.read_text()).get("library"))
    stamped.discard(None)
    # Joined rather than reduced to one name: a directory holding two libraries
    # is a thing worth seeing on the page, not a thing to pick a winner from.
    return " + ".join(sorted(stamped)) or None


def decay_heat_spread(results, material_id, case, steps, keys):
    """The nuclear-data sigma on the decay heat at each cooling step, in uW/g.

    Decay heat is a function of a whole inventory rather than of one nuclide, so
    the spread has to be taken over the ensemble's inventories rather than
    assembled from per-nuclide sigmas. A parent and its daughter move together
    under one resampled cross section, and adding their sigmas in quadrature
    would claim the independence the resampling exists precisely to avoid
    assuming.

    yani does that itself. ``get_decay_heat_uncertainty`` evaluates the heat once
    per replica and hands back the ensemble's sample standard deviation beside
    the unperturbed value, which is the same statistic this used to build by
    rebuilding a Material per replica out of ``get_uncertainty_inventories`` and
    taking ``np.std(..., ddof=1)`` over the totals. Both divide by N - 1, because
    the quantity wanted is how far the answer moves when the cross sections move
    within their covariance, which is the width of the ensemble rather than the
    error on where its centre sits.

    The getter works in watts and needs the material's volume, so everything is
    scaled here to the uW/g the measurement is published in.

    Returns the sigma on the total at each step and one dict of per-nuclide
    sigmas per step, both in uW/g, or ``None`` if the solve carried no ensemble.
    """
    totals = [results.get_decay_heat_uncertainty(material_id, step)
              for step in steps]

    # ``None`` for a solve run without ``data_uncertainty``, and a ``std_dev``
    # of ``None`` for one whose ensemble came back under two replicas, which is
    # what happens when nothing in the material had covariance to resample from:
    # cross sections converted with --no-covariance, or by a converter old
    # enough not to have written one, carry no covariance.arrow. A spread over
    # nothing is not zero and is not a number, and taking it anyway put a bare
    # NaN in the JSON, which is not valid JSON at all. Treat it as the absence
    # it is; `report_uncertainty` says so and names the fix.
    if any(estimate is None or estimate.std_dev is None for estimate in totals):
        return None

    scale = 1e6 / case.mass_g
    by_nuclide = []
    for step in steps:
        per_nuclide = results.get_decay_heat_uncertainty(
            material_id, step, by_nuclide=True)
        # Keyed over the same nuclides the breakdown carries, so a product the
        # ensemble never made lines up with a zero rather than dropping out of
        # the JSON on some steps and not others.
        by_nuclide.append({
            key: (per_nuclide[key].std_dev or 0.0) * scale if key in per_nuclide
            else 0.0
            for key in keys
        })
    return (np.array([estimate.std_dev * scale for estimate in totals]),
            by_nuclide)


def report_uncertainty(info, relative):
    """What the sigma covered, printed under the C/E it qualifies.

    A sigma is only readable next to what it left out. A zero on a nuclide whose
    evaluation states no MF=33 looks exactly like a zero on one that is well
    known, and the printed median would read as a whole-of-answer uncertainty
    when it covers the activation cross sections alone. Both come off the info
    yani files rather than being restated here, so this cannot drift from what
    the run actually did.
    """
    if relative is None:
        # The run asked for an uncertainty and got no ensemble back, which is
        # not the same as not asking. Say which, and say what to do about it.
        print("uncertainty: asked for, but no covariance was available to "
              "resample from.")
        print("             These cross sections carry no covariance.arrow. "
              "Reconvert them")
        print("             without --no-covariance to fill the +/- and %dCnuc "
              "columns:")
        print("             python convert_to_arrow.py --case <foil>")
        missing = info.get("no_covariance_data") or []
        if missing:
            print(f"             ({' '.join(sorted(missing))} reported no "
                  f"covariance data)")
        return

    lines = []
    if len(relative):
        lines.append(f"{np.median(relative) * 100:.1f}% median on the decay heat "
                     f"({relative.min() * 100:.1f}% to {relative.max() * 100:.1f}%)")

    samples, converged = info.get("samples"), info.get("converged")
    if samples is not None:
        settled = "converged" if converged else "hit the replica cap unconverged"
        lines.append(f"{samples} replicas, {settled}")

    perturbed = info.get("perturbed") or []
    missing = info.get("no_covariance_data") or []
    lines.append(f"{len(perturbed)} nuclides perturbed from MF=33 covariance")

    # Beside the count, the share of the production it actually reaches. A
    # library can state covariance for every isotope in the foil and still say
    # nothing about the one channel carrying the heat, and then the count reads
    # as coverage the sigma does not have. yani weights it by rate and by the
    # parent's own density, which this used to do by hand and got wrong: an edge
    # rate is per atom of its parent and production is not.
    fraction = info.get("rate_fraction_covered_total")
    if fraction is not None:
        lines.append(f"covering {fraction * 100:.1f}% of the production rate "
                     f"this irradiation drove"
                     + ("" if fraction > 0.9 else
                        ", so the sigma above is a spread over that part alone"))

    if missing:
        # Not a failure and not a warning to be silenced: it is the reason a
        # product can come back with a sigma of zero, which otherwise reads as
        # the most confident number in the table.
        lines.append(f"nothing published for {' '.join(sorted(missing))}, so "
                     f"reactions on them contribute no sigma")

    # A resampled rate that came out negative is truncated at zero rather than
    # used, which is right (a negative reaction rate is not a physical draw) and
    # is not free: the truncated draws are all on one side, so the ensemble mean
    # sits above the unperturbed value. Reported when it is common enough to
    # matter, on the same terms the report page prints it under the table.
    sampled, floored = info.get("rates_sampled") or 0, info.get("rates_floored") or 0
    if sampled and floored / sampled > 0.01:
        lines.append(f"{100.0 * floored / sampled:.1f}% of sampled rates went "
                     f"negative and were truncated at zero, biasing the mean up")

    clipped = info.get("matrices_clipped") or 0
    if clipped:
        worst = info.get("worst_relative_clip") or 0.0
        matrices = "matrix" if clipped == 1 else "matrices"
        lines.append(f"{clipped} covariance {matrices} repaired to positive "
                     f"semi-definite, worst by {worst * 100:.1f}%")

    not_perturbed = info.get("not_perturbed") or []
    if not_perturbed:
        lines.append("not propagated: " + ", ".join(not_perturbed))

    label = "uncertainty: "
    for index, line in enumerate(lines):
        print((label if index == 0 else " " * len(label)) + line)


def run(case, cross_sections, chain, uncertainty=None):
    """Specific decay heat [uW/g] after each cooling step, and its breakdown."""
    yani.cross_section_data = str(cross_sections)
    yani.transmutation_decay_data = DECAY_LIBRARY

    # Off, rather than pointed at a library. Nothing in any of the 73 FNS foils
    # fissions, so a yield source would never be read, and naming one anyway put
    # a library on the report's provenance page that contributed nothing to the
    # answer. `False` is yani's own setting for this and is not the same as
    # leaving it unset: unset falls back to a default library, whereas off means
    # a fission rate that did need yields is refused when the burnup matrix is
    # built, naming the nuclide, rather than quietly losing its fission products.
    yani.transmutation_fission_yields = False

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
            applied[subsection] = {"library": library_of(chain), "path": str(chain)}
        else:
            setattr(yani, setting, DECAY_LIBRARY)
            print(f"{subsection}: {DECAY_LIBRARY} ({warning})")
            applied[subsection] = {"library": DECAY_LIBRARY, "path": None}

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

    # What every part of the network was read from, recorded so the report can
    # state it instead of assuming it. No run here is on one library: the
    # neutron side moves with `--library`, the decay side never does, because no
    # neutron evaluation carries half-lives, decay modes or decay energies.
    # Naming both is the only way a reader can tell which of the two a
    # library-to-library difference belongs to.
    provenance = {
        "cross_sections": {"library": library_of(cross_sections, nuclides),
                           "path": str(cross_sections)},
        "reactions": applied.get("reactions"),
        "branching": applied.get("branching"),
        "decay_data": {"library": DECAY_LIBRARY, "path": None},
        "fission_yields": {"library": None, "path": None, "off": True},
    }

    # A chain left over from a different foil carries no reactions for these
    # nuclides and used to solve to an inventory of nothing, which surfaced only
    # as zero decay heat at the end of a full run. yani 0.13.0 refuses it before
    # the solve and names both the material and what the chain covers, so the
    # sidecar this file used to write and check is gone.

    # `transmute` hands back a TransmutationResults keyed by material id, not a
    # plain list. Its index 0 is the initial composition, so `step_materials` is
    # the per-step view the slice below counts from. The foil is built without an
    # explicit id, which files it under 0.
    #
    # `data_uncertainty` resamples the activation cross sections from the MF=33
    # covariance step 1 wrote, folds each draw against this spectrum and
    # re-solves. The solver itself is untouched and only its input moves, so the
    # mean inventories are what they would have been; omitting the argument reads
    # no covariance and costs nothing.
    results = material.transmute(schedule, data_uncertainty=uncertainty)
    material_id = material.id or 0
    states = results.step_materials(material_id)

    # One state per schedule step and no initial entry, so the irradiation
    # pulses come first and the measurement starts after them.
    heat, breakdown = [], []
    for state in states[len(case.irradiation):]:
        contributions = state.decay_heat(by_nuclide=True)
        heat.append(sum(contributions.values()) / case.mass_g * 1e6)
        breakdown.append(contributions)
    heat = np.array(heat)

    # The same cooling points, indexed the way the results are rather than the
    # way `states` is: `step_materials` drops the initial composition and the
    # results keep it at 0, so a cooling point j sits at step j + 1 in the
    # results once the irradiation pulses are counted off.
    cooling_steps = range(len(case.irradiation) + 1, results.num_steps + 1)
    sigma, by_nuclide_sigma, info = None, None, None
    if uncertainty is not None:
        info = results.data_uncertainty_info
        ensemble = decay_heat_spread(
            results, material_id, case, cooling_steps,
            keys=sorted({n for step in breakdown for n in step}))
        if ensemble is not None:
            sigma, by_nuclide_sigma = ensemble

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
        for parent, kinds in (results.get_reaction_rates(material_id, step)
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
    initial = results.get_material(material_id, 0)
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
    # The production routes into each product that carries heat, and the
    # flux-weighted isomeric branching, both derived by yani from the chain and
    # the rates the solve actually ran with. They were rebuilt in step 5 from
    # `edge_rates` plus a second load of the chain, which is a few hundred lines
    # of walking and one chance for the two chains to differ; yani 0.13.0
    # answers from the solve itself.
    #
    # Routes are asked for over the first irradiation pulse. Every pulse in
    # these schedules shares one spectrum, so the shares are the same for each,
    # and a cooldown drives no reactions at all.
    #
    # Only for products that carry heat. Walking every nuclide the chain can
    # reach costs time to answer a question about nuclides nobody will read a
    # row for; the floor is under step 5's own so the report can still choose.
    peak = {}
    for step in breakdown:
        total = sum(step.values()) or 1.0
        for nuclide, watts in step.items():
            peak[nuclide] = max(peak.get(nuclide, 0.0), watts / total)
    routes = {
        nuclide: results.get_production_routes(
            material_id, nuclide, 0,
            reaction_depth=ROUTE_REACTION_STEPS,
            decay_depth=ROUTE_DECAY_STEPS) or []
        for nuclide, share in sorted(peak.items())
        if share >= ROUTE_PRODUCT_FLOOR
    }
    branching = results.get_isomeric_branching(material_id, 0) or {}

    return (heat, breakdown, edge_rates, initial_atoms, sigma, by_nuclide_sigma,
            info, routes, branching, provenance)


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


def plot(case, calculated, breakdown, ratio, metrics, out_dir, sigma=None):
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

    # Both curves carry an uncertainty and the comparison is only readable when
    # both are drawn: a calculation half a decade from the measurement is a
    # different statement depending on whether its own band is 5% or 50% wide.
    if sigma is not None:
        handles.append(top.fill_between(
            times, calculated - sigma, calculated + sigma, color="#d62728",
            alpha=0.2, linewidth=0, zorder=1, label="cross-section covariance"))

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
    # Scaled to the two series being compared, the calculation and the
    # measurement, so that the comparison fills the panel. Short-lived
    # contributors fall away by orders of magnitude, and letting them set the
    # floor stretches the axis until the total and the measurement are one line.
    # They clip instead; the JSON carries every one of them at every step.
    edges = np.concatenate([calculated, measured, measured - uncertainty,
                            measured + uncertainty]
                           + ([calculated - sigma, calculated + sigma]
                              if sigma is not None else []))
    edges = edges[np.isfinite(edges) & (edges > 0)]
    if edges.size:
        top.set_ylim(bottom=edges.min() / 2.0, top=edges.max() * 2.0)
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
    if sigma is not None:
        relative = np.divide(sigma, calculated,
                             out=np.zeros_like(sigma), where=calculated > 0)
        bottom.fill_between(times, ratio * (1 - relative), ratio * (1 + relative),
                            color="#d62728", alpha=0.2, linewidth=0,
                            label="cross-section covariance")
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
                            "(default: data/<source>/chain-<case>, or the sibling "
                            "of --cross-sections when it ends with /neutron)")
    parser.add_argument("--output", type=pathlib.Path, default=None,
                        help="where the results go (default: results/<source>)")
    parser.add_argument("--list", action="store_true",
                        help="list the available foils and experiments, then exit")
    parser.add_argument("--source", default=None,
                        help="folder name to file this run's data and results under "
                             "(default: the library name, or the --endf-dir directory name)")
    parser.add_argument("--no-uncertainty", action="store_true",
                        help="skip the nuclear-data uncertainty, leaving the decay "
                             "heat as a bare value. The mean is identical either way")
    parser.add_argument("--samples", type=int, default=None,
                        help="fixed replica count for the uncertainty (default: let "
                             "yani add replicas until the sigmas settle). This bounds "
                             "cost or reproduces a specific run; it is not an "
                             "accuracy dial")
    parser.add_argument("--seed", type=int, default=42,
                        help="base seed for the resampling (default: 42). A nuclide's "
                             "perturbation is a function of (seed, replica, nuclide), "
                             "so the same seed reproduces a run whatever the replica "
                             "count")
    args = parser.parse_args()

    # Results and inputs live under the name of the nuclear data they came from,
    # so two libraries do not overwrite each other (see data_source).
    source = data_source.slug(args.library, args.endf_dir, args.source)
    if args.cross_sections is None:
        args.cross_sections = HERE / "data" / source / "neutron"
    if args.chain is None:
        args.chain = inferred_chain_path(args.cross_sections, args.case)
        if args.chain is None:
            args.chain = HERE / "data" / source / f"chain-{args.case}"
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

    uncertainty = None if args.no_uncertainty else yani.DataUncertainty(
        seed=args.seed, samples=args.samples)

    (calculated, breakdown, edge_rates, initial_atoms, sigma, by_nuclide_sigma,
     info, routes, branching, provenance) = run(
        case, args.cross_sections, args.chain, uncertainty)
    ratio, metrics = summarise(case, calculated)
    plot(case, calculated, breakdown, ratio, metrics, args.output, sigma)

    relative = (np.divide(sigma, calculated, out=np.zeros_like(sigma),
                          where=calculated > 0) if sigma is not None else None)
    header = f"\n{case.time_unit:>9}  {'measured':>10}  {'yani':>10}"
    header += "  +/-" if sigma is not None else ""
    print(header + f"  {'C/E':>6}")
    for index, (time, exp, calc, r) in enumerate(
            zip(case.times, case.measured, calculated, ratio)):
        line = f"{time:9.2f}  {exp:10.3e}  {calc:10.3e}"
        if relative is not None:
            line += f"  {relative[index] * 100:3.0f}%"
        print(line + f"  {r:6.3f}")
    print(f"\nmedian C/E {metrics['median_ratio']:.3f}, "
          f"mean deviation {metrics['mean_deviation_percent']:.1f}%, "
          f"measurement sigma {metrics['median_measurement_sigma_percent']:.1f}%")
    if info is not None:
        report_uncertainty(info, relative)

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
        # The nuclear-data uncertainty, absolute and in the same units as the
        # value beside it, so a reader never has to know which of the two a
        # percentage was taken against. `data_uncertainty_info` travels with it
        # because a sigma of zero means "no covariance published" as often as it
        # means "well known", and only the info says which.
        "yani_uncertainty_uW_per_g": None if sigma is None else sigma.tolist(),
        "by_nuclide_uncertainty_uW_per_g": by_nuclide_sigma,
        "data_uncertainty_info": info,
        # Derived by yani from the chain and the rates this solve ran with, so
        # the report binds them rather than deriving them from a chain it loads
        # separately and hopes is the same one.
        "production_routes": routes,
        "isomeric_branching": branching,
        # Which library each part of the network came from. Filed with the
        # result rather than worked out later from the paths, because the paths
        # can be repointed and this run cannot be run again from its own JSON.
        "nuclear_data": provenance,
    }, indent=2))


if __name__ == "__main__":
    main()
