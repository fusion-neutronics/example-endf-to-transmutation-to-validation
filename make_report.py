#!/usr/bin/env python3
"""Step 5: bind the sweep results into a decay-heat validation report PDF.

Reads the per-foil JSON that ``run_transmutation.py`` and ``sweep_fns.py``
already write, one sweep per library, and lays them out in the page order the
UKAEA decay-heat validation reports use:

    make_report.py --case W --libraries tendl-2025 tendl-2017

    cover page               the foil, named
    C/E table                one row per cooling point, one E/C column per library
    dominant products        what carries the heat, with half-lives
    production pathways      which reaction on which parent makes each product
    figure page              heat curves and % contributions, one row per library

Nothing is recomputed: the inventory solve happened in step 2 and its full
per-nuclide breakdown is in the JSON, so a report is cheap to regenerate and
cannot disagree with the run it came from.

Three columns the published reports carry are NOT reproduced here, because the
saved results do not contain what they need. They are left out rather than
filled with a plausible number:

* An uncertainty on the calculated value (the reports print "+/- 6%"). That is
  cross-section covariance propagated through the inventory, which yani does
  not currently produce, so the calculated column here is a bare value.
* "Path %", the fraction of a product's inventory that arrives down each
  pathway. That needs per-reaction rates, sigma(g) folded with the group flux
  for every edge, and the results JSON keeps only the inventory the rates
  produced. The pathways themselves are real, taken from the library's own
  topology; only their relative weights are missing.
* Decay feeds beyond one neutron reaction. The chain exposes half-lives and the
  reaction topology but not decay modes, so ``W186(n,2n)W185m`` is derivable and
  the ``(IT)W185`` step that follows it is not. Isomers are marked as feeding
  their ground state, which is the one decay edge a name settles on its own.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import tempfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent

# A4 portrait, which is what the published reports are set on.
PAGE = (8.27, 11.69)

# The rule under the running header.
RULE = "#e8a33d"

# Monospace so the tables column up without measuring glyphs.
MONO = {"family": "monospace", "fontsize": 7.4}

# A product is named on the pathway and product pages when it carries at least
# this share of the total heat at any one cooling point. A fixed count cannot
# work across cases: the 5-minute foils run on a handful of short-lived isomers
# and the long-cooling ones on a different, smaller set, and the published
# reports plot 10 and 6 nuclides respectively for exactly that reason.
#
# Note this is a peak share, not a share of the heat summed over the cooling
# steps. Summed heat is dominated by the early high-heat points, so ranking on
# it buries any product that peaks late: on W/2000exp_5min it puts Ta184 tenth
# despite Ta184 carrying 18% of the heat at one point.
SHARE_FLOOR = 0.002

# The cover page names the foil. Element symbols only; the benchmark's four
# alloy foils (Inc600, NiCr, SS304, SS316) are already their own names.
# The decay library the reaction data is paired with. TENDL publishes no decay
# sublibrary, so this is the same default run_transmutation.py holds, and the
# same reason: half-lives have to come from somewhere.
DECAY_LIBRARY = "endf-b8.1"

# Where convert_to_arrow.py's cache puts a converted decay sublibrary.
CACHE = pathlib.Path.home() / ".cache" / "yamc"


def element_name(case):
    """The foil's name for the cover page.

    The benchmark's four alloy foils (Inc600, NiCr, SS304, SS316) are already
    their own names and are returned unchanged.
    """
    import yani.data

    return yani.data.element_names().get(case, case).title()


def foil_isotopes(case, experiment):
    """The foil's own natural isotopes, which are the only pathway starts.

    Read from the benchmark deck rather than inferred from the foil name, since
    a quarter of the alloy foils would give the wrong answer either way.
    """
    import fns_case
    import yani.data

    natural = yani.data.element_nuclides()
    composition = fns_case.load(case, experiment).composition
    return [nuclide for element in sorted(composition)
            for nuclide in natural.get(element, [])]


def read_chain(chain, decay=DECAY_LIBRARY, cache=CACHE):
    """The chain, with a decay subsection borrowed if it has none of its own.

    convert_to_arrow.py writes a chain as ``reactions/`` and ``branching/``, and
    leaves decay data to the converted sublibrary in the cache, which is why
    run_transmutation.py points two settings at two places. TransmutationChain
    wants one directory holding all of it, so when the chain has no ``decay/``
    of its own it is composed here: the chain's own subsections, plus a decay
    subsection linked in from the cache.

    Composing it rather than falling back to the cache's own chain is the whole
    point. The cached chain's topology sends every (n,2n) to a ground state,
    because the isomeric branching is energy dependent and lives in the
    ``branching/`` subsection the converter wrote. Without it W186(n,2n)W185m
    does not exist, and on a tungsten foil that is 98% of the decay heat.
    """
    if chain is None:
        return None, None
    try:
        import yani

        if (chain / "decay").is_dir():
            return yani.TransmutationChain(str(chain)), None

        borrowed = cache / f"{decay}-transmutation-decay.arrow"
        if not borrowed.is_dir():
            print(f"chain: {chain} has no decay/ and {borrowed} is not there "
                  f"either, so half-lives and pathways are unavailable.\n"
                  f"       Convert the decay sublibrary first, or pass --chain a "
                  f"directory that has its own decay/.")
            return None, None

        # Held open for as long as the chain is read, and removed after.
        composed = tempfile.TemporaryDirectory(prefix="yani-report-chain-")
        root = pathlib.Path(composed.name)
        subsections = {"decay": borrowed}
        for name in ("reactions", "branching", "fission_yields"):
            if (chain / name).is_dir():
                subsections[name] = (chain / name).resolve()
        for name, target in subsections.items():
            (root / name).symlink_to(target, target_is_directory=True)
        (root / "manifest.json").write_text(json.dumps({
            "format_version": 2,
            "library": f"{chain.name}+{decay}",
            "subsections": {name: {"path": name} for name in subsections},
        }))
        loaded = yani.TransmutationChain(str(root))
        print(f"chain: {chain} with decay from {borrowed.name} "
              f"({', '.join(sorted(subsections))})")
        return loaded, composed
    except Exception as error:  # noqa: BLE001 - both pages that use it are optional
        print(f"chain: unreadable at {chain} ({error}). Half-lives and pathways "
              f"are left out.")
        return None, None


def topology_source(chain, seeds, depth, case):
    """Reaction topology as target -> [(parent, reaction)], scoped to the foil.

    Two scopes are applied, and both matter.

    The chain covers every nuclide the library knows, so asking it what makes
    W187 answers with Os190(n,a) and Ir192(n,npa) as readily as W186(n,gamma).
    Those are real edges and irrelevant here: nothing in a tungsten foil is
    osmium. So the walk starts at the foil's own isotopes.

    Walking from there to saturation does not fix it, because successive capture
    does reach osmium and iridium from tungsten eventually. It takes a fluence
    the FNS foils never see: these irradiations are 5 minutes to 7 hours, where
    the inventory is first and second order in the reaction rate and nothing
    else registers. ``depth`` is the cap that keeps the page to that, and it
    counts reaction steps, so decay does not spend it.
    """
    if chain is None:
        return {}, ("No chain was given, so the topology is unavailable. Pass "
                    "--chain the\ndirectory convert_to_arrow.py wrote for this "
                    "foil to fill this page.")
    forward = chain.reactions

    # A chain is scoped to the foil it was converted for, so one built for a
    # different foil loads perfectly and carries no reaction for anything in
    # this one. That is worth saying, because it is not the same problem as
    # having no chain and it has a different fix. run_transmutation.py refuses
    # outright on this; here it costs one page, so it reports instead.
    covered = [seed for seed in seeds if forward.get(seed)]
    if not covered:
        return {}, ("The chain given carries no reaction for any isotope in this "
                    "foil, so it was\nbuilt for a different one. Rebuild it for "
                    "this foil:\n"
                    "  python convert_to_arrow.py --case " + case +
                    " --chain <chain>")

    # Breadth first from the foil's isotopes, so a nuclide's depth is the fewest
    # reaction steps that reach it, and an edge is kept only if its parent sits
    # inside the cap.
    reached = {seed: 0 for seed in covered}
    frontier = list(reached)
    for step in range(1, depth + 1):
        following = []
        for parent in frontier:
            for _reaction, target, _ratio in forward.get(parent, []):
                if target and target not in reached:
                    reached[target] = step
                    following.append(target)
        frontier = following

    backward = {}
    for parent, level in reached.items():
        if level >= depth:
            continue
        for reaction, target, _ratio in forward.get(parent, []):
            if target:
                backward.setdefault(target, []).append((parent, reaction))
    return backward, None


def edge_order(abundances):
    """Sort key putting a product's parents in descending natural abundance.

    Not a contribution ranking, which would need the reaction rates, but not
    arbitrary either: a route off a parent that is 28% of the foil is a better
    first line than one off a parent that is 0.12%, and W180 being 0.12% of
    natural tungsten is exactly why W179 is a smaller product than W185.
    """
    def key(edge):
        parent, reaction = edge
        return -abundances.get(parent, 0.0), parent, reaction

    return key


def format_half_life(seconds):
    """A half-life the way the decay tables print it."""
    if seconds is None:
        return "-"
    if seconds <= 0:
        return "stable"
    # Thresholds are deliberately not the unit sizes: 23.83 h is how the decay
    # tables print W187, and rolling it into "1d" at 86400 s loses the figure
    # the table exists to show.
    for size, unit, floor in ((86400 * 365.25, "y", 86400 * 365.25),
                              (86400, "d", 86400 * 3),
                              (3600, "h", 3600),
                              (60, "m", 60),
                              (1, "s", 0)):
        if seconds >= floor:
            return f"{seconds / size:.4g}{unit}"
    return f"{seconds:.3g}s"


def calculated_of(result):
    """The calculated heat array, under whichever key wrote it.

    ``run_transmutation.py`` writes ``yani_uW_per_g``. Results filed before the
    library was renamed carry the same array under ``yats_uW_per_g``, and those
    sweeps are expensive enough to be worth reading rather than repeating.
    """
    for key in ("yani_uW_per_g", "yats_uW_per_g"):
        if key in result:
            return np.array(result[key], dtype=float)
    raise SystemExit(
        f"{result.get('case')}/{result.get('experiment')}: no calculated heat in "
        f"the result. Expected a yani_uW_per_g key; found {sorted(result)}."
    )


def load_results(case, experiment, libraries, results_root):
    """One result per library for this foil, in the order the libraries were given.

    A library with no result for the foil is reported and dropped rather than
    failing the report: a converter that could not do W has not thereby made the
    other libraries unreportable.
    """
    found, missing = [], []
    for library in libraries:
        for root in (results_root / library / "sweep", results_root / library):
            path = root / f"fns_{case}_{experiment}.json"
            if path.is_file():
                found.append((library, json.loads(path.read_text())))
                break
        else:
            missing.append(library)
    if missing:
        print(f"no {case}/{experiment} result for: {', '.join(missing)}")
    if not found:
        raise SystemExit(
            f"no {case}/{experiment} results under {results_root} for any of "
            f"{', '.join(libraries)}.\n"
            f"Run: python run_transmutation.py --case {case} --source <source> "
            f"--output {results_root / libraries[0]}"
        )
    return found


def score(result):
    """The two figures of merit the published tables carry, per library.

    ``mean % diff. from E`` is the mean of |C/E - 1|, and ``mean chi^2`` divides
    by the measurement's own sigma, so it says whether a deviation is larger
    than the experiment can resolve. Both are the report's definitions: neither
    carries an uncertainty on the calculated value, because there is not one.
    """
    measured = np.array(result["measured_uW_per_g"], dtype=float)
    sigma = np.array(result["measured_uncertainty"], dtype=float)
    calculated = calculated_of(result)
    good = measured > 0
    ratio = np.divide(calculated, measured, out=np.full_like(calculated, np.nan),
                      where=good)
    chi2 = np.divide((calculated - measured) ** 2, sigma ** 2,
                     out=np.full_like(calculated, np.nan), where=sigma > 0)
    return {
        "measured": measured,
        "sigma": sigma,
        "calculated": calculated,
        "ratio": ratio,
        "mean_percent_diff": float(np.nanmean(np.abs(ratio - 1.0)) * 100.0),
        "mean_chi2": float(np.nanmean(chi2)),
    }


def leading_products(result, floor=SHARE_FLOOR):
    """Products carrying at least ``floor`` of the total heat at some one point.

    Ordered by that peak share, so the ordering is "how big does this ever get"
    rather than "how big is this while the total is largest".
    """
    steps = result["by_nuclide_uW_per_g"]
    totals = [sum(step.values()) for step in steps]
    peak, when = {}, {}
    for index, (step, total) in enumerate(zip(steps, totals)):
        if total <= 0:
            continue
        for nuclide, heat in step.items():
            share = heat / total
            if share > peak.get(nuclide, 0.0):
                peak[nuclide], when[nuclide] = share, index
    keep = [n for n, share in peak.items() if share >= floor]
    keep.sort(key=lambda n: -peak[n])
    return [(n, peak[n], when[n]) for n in keep]


def series_of(result, nuclide):
    """One product's specific heat through the cooling steps."""
    return np.array([step.get(nuclide, 0.0) for step in result["by_nuclide_uW_per_g"]],
                    dtype=float)


def page(pdf, title, subtitle, number, total):
    """A blank page carrying the running header, the rule and the footer."""
    figure = plt.figure(figsize=PAGE)
    figure.text(0.94, 0.962, title, ha="right", va="bottom", fontsize=9)
    figure.text(0.94, 0.944, subtitle, ha="right", va="bottom", fontsize=8,
                style="italic")
    figure.add_artist(plt.Line2D([0.30, 0.94], [0.938, 0.938], color=RULE,
                                 linewidth=1.1, transform=figure.transFigure))
    figure.text(0.5, 0.028, f"Page {number} of {total}", ha="center", fontsize=8)
    return figure


def block(figure, lines, top, left=0.08, pitch=0.0132, **style):
    """Fixed-pitch text rows, returning the height left below them.

    The tables are drawn as text rather than assembled with ax.table, because
    twenty rows of numbers want a monospace column and a known line pitch, not
    cell autoscaling.
    """
    font = dict(MONO)
    font.update(style)
    for offset, line in enumerate(lines):
        figure.text(left, top - offset * pitch, line, va="top", **font)
    return top - len(lines) * pitch


def cover_page(pdf, case, results, title, subtitle, number, total):
    """The foil, named, as the published sections open."""
    figure = page(pdf, title, subtitle, number, total)
    figure.text(0.5, 0.62, element_name(case), ha="center", fontsize=26)
    _library, first = results[0]
    figure.text(0.5, 0.56, f"FNS decay heat validation, {first['experiment']}",
                ha="center", fontsize=11)
    figure.text(0.5, 0.53, f"libraries: {', '.join(lib for lib, _ in results)}",
                ha="center", fontsize=9, color="#555555")
    pdf.savefig(figure)
    plt.close(figure)


def table_page(pdf, case, results, half_lives, title, subtitle, number, total):
    """The C/E table: measured, the primary library's value, then E/C per library.

    E/C rather than C/E, to read against the published tables. The repo's JSON
    stores C/E, so this is the reciprocal of what step 2 printed.
    """
    figure = page(pdf, title, subtitle, number, total)
    primary_name, primary = results[0]
    scores = [(name, score(result)) for name, result in results]
    unit = primary["time_unit"]
    times = np.array(primary["times"], dtype=float)

    head = (f"{'Times':>8} {'FNS Exp.':>22} {primary_name[:14]:>14}"
            + "".join(f"{name[:12]:>13}" for name, _ in scores))
    units = (f"{unit[:8]:>8} {'uW/g':>22} {'uW/g':>14}"
             + "".join(f"{'E/C':>13}" for _ in scores))
    lines = [head, units, "-" * len(head)]
    for index, moment in enumerate(times):
        measured = primary["measured_uW_per_g"][index]
        sigma = primary["measured_uncertainty"][index]
        percent = 100.0 * sigma / measured if measured else float("nan")
        row = (f"{moment:8.2f} {measured:11.2E} +/- {percent:3.0f}% "
               f"{scores[0][1]['calculated'][index]:14.2E}")
        for _name, values in scores:
            ratio = values["ratio"][index]
            row += f"{1.0 / ratio:13.2f}" if ratio and np.isfinite(ratio) else f"{'-':>13}"
        lines.append(row)
    lines.append("-" * len(head))
    lines.append(f"{'mean % diff. from E':>46}"
                 + "".join(f"{v['mean_percent_diff']:13.0f}" for _, v in scores))
    lines.append(f"{'mean chi^2':>46}"
                 + "".join(f"{v['mean_chi2']:13.2f}" for _, v in scores))

    figure.text(0.08, 0.905, f"{element_name(case)}, {primary['experiment']}",
                fontsize=11)
    bottom = block(figure, lines, 0.878)

    # The products, with their half-lives, so the table above is readable as a
    # decay curve rather than as twenty numbers.
    products = leading_products(primary)
    figure.text(0.08, bottom - 0.022,
                f"{primary_name} dominant products, by peak share of the total",
                fontsize=10)
    rows = [f"{'Product':<10}{'T1/2':>10}{'peak share':>12}{'at time':>12}"
            f"{'first point':>14}",
            "-" * 58]
    first_total = sum(primary["by_nuclide_uW_per_g"][0].values())
    for nuclide, share, index in products:
        at_first = primary["by_nuclide_uW_per_g"][0].get(nuclide, 0.0)
        rows.append(
            f"{nuclide:<10}{format_half_life(half_lives.get(nuclide)):>10}"
            f"{share * 100:11.1f}%{times[index]:11.2f} {unit[:1]}"
            f"{(at_first / first_total * 100 if first_total else 0.0):13.1f}%"
        )
    block(figure, rows, bottom - 0.040)
    pdf.savefig(figure)
    plt.close(figure)


def pathway_page(pdf, case, results, topology, note, half_lives, depth, abundances,
                 title, subtitle, number, total):
    """Which reaction on which parent makes each product that carries heat."""
    figure = page(pdf, title, subtitle, number, total)
    primary_name, primary = results[0]
    figure.text(0.08, 0.905,
                f"{element_name(case)}, {primary_name} production pathways",
                fontsize=11)
    seen = {nuclide for nuclide, _share, _at in leading_products(primary)}

    if not topology:
        figure.text(0.08, 0.86, note or "No topology was available.",
                    va="top", fontsize=8.6, family="monospace")
        pdf.savefig(figure)
        plt.close(figure)
        return

    rows = [f"{'Product':<10}{'T1/2':>10}  {'Pathway':<38}{'peak share':>12}",
            "-" * 72]
    for nuclide, share, _index in leading_products(primary):
        feeds = [f"{parent}{reaction}{nuclide}"
                 for parent, reaction in sorted(topology.get(nuclide, []),
                                                key=edge_order(abundances))]
        # An isomer feeding its own ground state is the one decay edge a name
        # settles without a decay mode, and it is worth showing because it is
        # often the larger half of how the ground state is made: the report has
        # W186(n,2n)W185m(IT)W185 carrying more of W185 than the direct
        # W186(n,2n)W185 does. It belongs to the ground state, not to the
        # isomer, so it is added here rather than under the isomer's own row.
        for isomer in (f"{nuclide}_m1", f"{nuclide}_m2"):
            if isomer in topology or isomer in seen:
                feeds.append(f"{isomer}(IT){nuclide}")
        if not feeds:
            feeds = ["no pathway within the depth searched"]
        rows.append(f"{nuclide:<10}{format_half_life(half_lives.get(nuclide)):>10}  "
                    f"{feeds[0]:<38}{share * 100:11.1f}%")
        rows.extend(f"{'':<22}{feed:<38}" for feed in feeds[1:])
    bottom = block(figure, rows, 0.878)
    figure.text(0.08, bottom - 0.020,
                f"Pathways are the library's own topology, walked {depth} reaction "
                f"step{'s' if depth != 1 else ''} from the foil's\n"
                "own isotopes, in descending parent abundance. That is not a ranking:\n"
                "the relative weights (the reports' \"Path %\") need per-reaction rates,\n"
                "and the saved inventory carries only what those rates produced.\n"
                "\"peak share\" is the product's own largest share of the total heat.",
                va="top", fontsize=7.6, color="#555555")
    pdf.savefig(figure)
    plt.close(figure)


def figure_page(pdf, case, results, title, subtitle, number, total):
    """Heat curves left, % contributions right, one row per library."""
    figure = page(pdf, title, subtitle, number, total)
    rows = len(results)
    # hspace is a fraction of subplot height and subplot height goes as 1/rows,
    # so a fixed hspace leaves a gap that grows as rows are removed. Scaling it
    # with the row count keeps the gap the same on a 2-library report as on a 4.
    # bottom clears the footer, which the axis label used to be drawn over.
    grid = figure.add_gridspec(rows, 2, left=0.09, right=0.90, top=0.90,
                               bottom=0.085, hspace=0.17 * rows, wspace=0.58)

    for row, (library, result) in enumerate(results):
        values = score(result)
        times = np.array(result["times"], dtype=float)
        unit = result["time_unit"]
        products = leading_products(result)
        # tab10 repeats after ten and these pages routinely name more than ten
        # products, which put two different nuclides on the same colour.
        palette = plt.get_cmap("tab20" if len(products) > 10 else "tab10")
        colours = [palette(i % palette.N) for i in range(len(products))]

        heat = figure.add_subplot(grid[row, 0])
        heat.errorbar(times, values["measured"], yerr=values["sigma"], fmt="^",
                      color="#808080", markersize=3, linewidth=0.7,
                      label="FNS Experiment", zorder=3)
        for index, (nuclide, _share, _at) in enumerate(products):
            heat.plot(times, series_of(result, nuclide), "--", linewidth=0.7,
                      color=colours[index], label=nuclide)
        heat.plot(times, values["calculated"], "-", color="black", linewidth=1.2,
                  label="Total", zorder=2)
        heat.set_xscale("log")
        heat.set_yscale("log")
        positive = values["calculated"][values["calculated"] > 0]
        if positive.size:
            heat.set_ylim(positive.min() / 1e3, positive.max() * 5)
        heat.set_xlabel(f"Time after irradiation [{unit}]", fontsize=6.5)
        heat.set_ylabel("Heat Output [uW/g]", fontsize=6.5)
        heat.set_title(f"FNS {result['experiment']} - {case} - {library}", fontsize=7.5)
        heat.tick_params(labelsize=6)
        heat.legend(fontsize=5.0, loc="upper left", bbox_to_anchor=(1.01, 1.0),
                    frameon=False, handlelength=1.5, borderpad=0.2,
                    labelspacing=0.28)
        heat.grid(alpha=0.25, linewidth=0.4)

        share = figure.add_subplot(grid[row, 1])
        totals = np.array([sum(step.values()) for step in result["by_nuclide_uW_per_g"]])
        for index, (nuclide, _peak, _at) in enumerate(products):
            percent = np.divide(series_of(result, nuclide) * 100.0, totals,
                                out=np.zeros_like(totals), where=totals > 0)
            share.plot(times, percent, "--", linewidth=0.7,
                       color=colours[index], label=nuclide)
        share.set_xscale("log")
        share.set_ylim(0, 100)
        share.set_xlabel(f"Time after irradiation [{unit}]", fontsize=6.5)
        share.set_ylabel("%decay heat contributions", fontsize=6.5)
        share.set_title(f"FNS {result['experiment']} - {case} - {library}", fontsize=7.5)
        share.tick_params(labelsize=6)
        share.legend(fontsize=5.0, loc="upper left", bbox_to_anchor=(1.01, 1.0),
                     frameon=False, handlelength=1.5, borderpad=0.2,
                     labelspacing=0.28)
        share.grid(alpha=0.25, linewidth=0.4)

    pdf.savefig(figure)
    plt.close(figure)


def build(case, experiment, libraries, results_root, chain, decay, depth, out,
          title, subtitle):
    """Write the report and return where it went."""
    import yani.data

    results = load_results(case, experiment, libraries, results_root)
    loaded, composed = read_chain(chain, decay)
    half_lives = dict(loaded.half_lives) if loaded else {}
    topology, note = topology_source(loaded, foil_isotopes(case, experiment), depth,
                                     case)
    if note:
        print(f"pathways: {note.splitlines()[0]}")
    abundances = yani.data.natural_abundance()

    # Four pages, so the footer can count them before any is drawn.
    total = 4
    out.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(out) as pdf:
        cover_page(pdf, case, results, title, subtitle, 1, total)
        table_page(pdf, case, results, half_lives, title, subtitle, 2, total)
        pathway_page(pdf, case, results, topology, note, half_lives, depth,
                     abundances, title, subtitle, 3, total)
        figure_page(pdf, case, results, title, subtitle, 4, total)
        pdf.infodict()["Title"] = f"{title}: {element_name(case)}"
        pdf.infodict()["Subject"] = subtitle

    if composed is not None:
        composed.cleanup()

    for library, result in results:
        values = score(result)
        print(f"{library:>14}: mean % diff {values['mean_percent_diff']:6.1f}   "
              f"mean chi^2 {values['mean_chi2']:8.2f}   "
              f"{len(leading_products(result))} products above "
              f"{SHARE_FLOOR * 100:g}% peak share")
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--case", required=True,
                        help="foil to report, as named in the benchmark (W, Fe, SS316)")
    parser.add_argument("--experiment", default="2000exp_5min",
                        help="which experiment's results to read "
                             "(default: 2000exp_5min)")
    parser.add_argument("--libraries", nargs="+", required=True,
                        help="folders under --results, one per library. The first "
                             "is the primary one, whose absolute values the table "
                             "carries")
    parser.add_argument("--results", type=pathlib.Path, default=HERE / "results",
                        help="root the sweeps were filed under (default: results/)")
    parser.add_argument("--chain", type=pathlib.Path, default=None,
                        help="chain directory convert_to_arrow.py wrote, for "
                             "half-lives and the pathway page. Its branching/ is "
                             "what puts isomers on the page, so prefer the "
                             "library's own over a generic chain. Optional: "
                             "without it those columns are left out")
    parser.add_argument("--decay", default=DECAY_LIBRARY,
                        help="decay sublibrary supplying half-lives when --chain has "
                             f"no decay/ of its own (default: {DECAY_LIBRARY})")
    parser.add_argument("--pathway-depth", type=int, default=2,
                        help="how many reaction steps from the foil's own isotopes "
                             "the pathway page walks. These irradiations are short "
                             "enough that the inventory is first and second order, "
                             "so 2 is the honest depth (default: 2)")
    parser.add_argument("--output", type=pathlib.Path, default=None,
                        help="PDF to write (default: "
                             "results/report_<case>_<experiment>.pdf)")
    parser.add_argument("--title", default="yani decay heat validation",
                        help="running header, top line")
    parser.add_argument("--subtitle", default="FNS benchmark, JAEA / IAEA CoNDERC",
                        help="running header, second line")
    args = parser.parse_args()

    out = args.output or (args.results / f"report_{args.case}_{args.experiment}.pdf")
    path = build(args.case, args.experiment, args.libraries, args.results,
                 args.chain, args.decay, args.pathway_depth, out, args.title,
                 args.subtitle)
    print(f"\nreport: {path}")


if __name__ == "__main__":
    main()
