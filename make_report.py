#!/usr/bin/env python3
"""Step 5: bind the sweep results into a decay-heat validation report PDF.

Reads the per-foil JSON that ``run_transmutation.py`` and ``sweep_fns.py``
already write, one sweep per library, and lays them out in the page order a
decay-heat validation report conventionally uses:

    make_report.py                 # one document per foil, all of its campaigns
    make_report.py --case W        # just this one

    cover page               the foil, named, and one line per campaign
    library comparison       every library's total against the measurement, one
                             panel per campaign, each curve carrying its own
                             cross-section uncertainty band
    C/E table                one row per cooling point, one E/C column per
                             library, the nuclide E/C analysis under it saying
                             which product carries the disagreement, and that
                             library's heat curve under that
    production pathways      the reaction and decay steps that make each
                             product, over as many pages as they need
    figure page              heat curves and % contributions, one row per library

The last three repeat per campaign; the comparison page carries all of them at
once, so the campaigns can be read against each other. One foil is one
document, and that is the point of binding it this way. A foil
measured more than once is still one subject, and the spread between its
campaigns is a result about the data rather than about any one measurement:
iron reads 6% high against ``2000exp_5min`` and 7% low against
``1996exp_5min``, which is a statement neither report makes on its own.
Tungsten's three run 122%, 65% and 20% out.

Which foils, which campaigns, which libraries and which chain are all worked out
from what is on disk, so the common case takes no arguments at all. Each can be
given instead, and ``--libraries`` is also how to choose which library is the
primary one.

The same tables are written beside the PDF in a form something else can read: a
JSON carrying all of them, uncapped, and a CSV of the C/E table.

Nothing is recomputed: the inventory solve happened in step 2 and its full
per-nuclide breakdown is in the JSON, so a report is cheap to regenerate and
cannot disagree with the run it came from.

Three columns mean something the heading does not say on its own:

* An uncertainty on the calculated value, and with it ``%dC_nuc`` in the
  nuclide analysis. Both are cross-section covariance carried through the
  inventory. yani reads the MF=33 covariance step 1 writes beside each nuclide,
  resamples the activation cross sections from it, folds each draw against the
  foil's own spectrum and re-solves, so the spread over that ensemble is the
  number. Step 2 writes it into its JSON and both columns are read from there.

  The total is taken over whole inventories rather than assembled from
  per-nuclide sigmas. Decay heat is a function of an inventory, and a parent and
  its daughter move together under one resampled cross section, so quadrature
  over nuclides would claim an independence that the resampling exists to avoid
  assuming. ``%dC_nuc`` is the narrower thing: one product's own spread, which
  is what says whether the row beside it is inside what its cross sections can
  account for.

  A sigma of zero is not a claim of certainty. An evaluation carrying no MF=33
  is perturbed by nothing, so its products come back exact; step 2 files
  ``data_uncertainty_info`` alongside, which names those nuclides, and a run
  made without the covariance leaves both columns empty rather than at zero.

* "Path %", the fraction of a product's inventory arriving down each route.
  ``get_production_routes`` answers from the chain and the rates the solve ran
  with, weighted by what the foil's own nuclides drove, and step 2 files the
  answer. A route is followed past its neutron reaction, so
  ``W186(n,2n)W185_m1(IT)W185`` is read rather than guessed at.

* The flux-weighted isomeric branching, which says which evaluated quantity a
  disagreement belongs to. It cannot be read off the chain: the chain file's own
  branching for the dominant tungsten channel is a placeholder,

      W186 (n,2n) -> W185      branching 1.000000
      W186 (n,2n) -> W185_m1   branching 0.000000

  which the energy-dependent overlay in ``branching/`` replaces at solve time.
  It is not zero in any meaningful sense: with the overlay configured W185m
  carries 98% of the decay heat at the first cooling point, and without it 10%.
  The split the solve actually used is printed rather than inferred.

A result carrying no sigma leaves the two uncertainty columns empty rather than
at zero, and one carrying no routes leaves that page for a rerun of step 2.
What this script will not do is invent either from what is there.

One difference from the published pages is deliberate rather than a gap. They
write an isomer as ``W185m`` and a second metastable state as ``Ta182n``; every
table here writes ``W185_m1`` and ``Ta182_m2``, which is what yani calls them.
Keeping the two spellings the same means a name on the page is the key in the
result JSON and in the chain, so a row can be traced back to the data it came
from without a translation table in between.
"""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import tempfile
import textwrap

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent

# A4 portrait, which is what the published reports are set on.
PAGE = (8.27, 11.69)

# Panel geometry on the figure page. The width left between the two columns is
# for the legends, which sit outside the axes and routinely name more than ten
# products. The aspect is height over width, a little wider than square, which
# is what a decay curve over three decades of time wants.
PANEL_WSPACE = 0.58
PANEL_ASPECT = 0.78

# How tall a panel may grow when the page has height to spare. Two panels
# abreast fix the width, so a two-library page laid out at PANEL_ASPECT fills
# barely half the sheet; the rows grow into the slack up to this before any of
# it is left as margin. Past about here a panel starts to read as a column
# rather than a plot, which is the shape the fixed layout was drawn to avoid.
PANEL_ASPECT_MAX = 1.15

# How far past the data a heat panel's y axis runs, as a factor on each end.
# 2.0 is about a third of a decade of air above and below, enough that the top
# curve is not jammed against the frame and little enough that it is still the
# data setting the scale.
PANEL_HEADROOM = 2.0

# How solid the +/- 1 sigma band around a calculated curve is drawn. Light
# enough that four of them can overlap on the comparison panel and the curves
# stay readable through them, dark enough that a band is visibly a band and not
# a rendering artefact. The bands routinely do overlap: that they do is the
# point, since two libraries whose bands intersect do not disagree.
BAND_ALPHA = 0.18

# The share of the driven production rate a library's covariance has to span
# before its sigma is left to speak for itself. Below this the page says so
# under the table, because the number is then a spread over part of the answer
# and reads as a spread over all of it. Set just under 1 rather than at some
# round fraction: anything that leaves a tenth of the production unperturbed is
# worth a sentence, and the case this exists for spans 5%.
COVERAGE_FLOOR = 0.9

# Room below a panel's axes box for its tick labels and x label, as a fraction
# of the page. matplotlib places both outside the box, so anything laid out
# against the box's own bottom edge sits on top of them.
XLABEL_CLEAR = 0.042

# The least height, as a fraction of the page, worth giving the heat curve that
# sits under the C/E table. What the tables leave is not fixed: a foil with 20
# cooling points and 13 products leaves less than one with 10 points and 5. Any
# less than this and the decay spans three decades in under an inch, which is
# unreadable and is already drawn properly on the figure page, so the space is
# left as margin instead.
PANEL_MIN_HEIGHT = 0.17

# The rule under the running header.
RULE = "#e8a33d"

# Monospace so the tables column up without measuring glyphs.
MONO = {"family": "monospace", "fontsize": 7.4}

# Space kept clear at the right edge of a right-aligned table column, as a
# fraction of the page width, so its figures do not run into the next column.
GUTTER = 0.009

# Ink for the booktabs rules, and for the hairlines that separate the row
# groups inside a table. The hairline is light enough to read as a guide for
# the eye rather than as a rule dividing one table into several.
RULE_INK = "#333333"
FAINT_INK = "#bfbfbf"

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

# The share of the calculation a product needs to earn a line in the nuclide
# E/C analysis, which is a stricter question than "does this product carry
# heat". That table exists to say where a disagreement lives, and a product
# that is 0.3% of the calculation cannot account for a C/E of 1.9 whatever its
# cross section is doing. The published page prints one row for tungsten's 5
# minute irradiation; this prints the six that could each move the answer, and
# leaves the long tail to the pathway page and the figure legends, which is
# also what keeps room on the page for the curve underneath.
NUCLIDE_EC_FLOOR = 0.05

# Routes shown per product. Products reached by many small channels can carry
# dozens once decay steps are included, and a page of them buries the one that
# matters. Whatever this drops is counted and said on the page.
ROUTES_SHOWN = 6

# The share of a product's production below which its route is not printed.
# A cap alone is not enough: the cap fills its six lines whether or not there
# are six routes worth reading, so a product made almost entirely down one
# channel still gets five lines of 0.0% under it, and a page of those reads as
# if the routes were unranked. The published pages carry nothing under about
# 2%, which this is well under; it is set to drop what rounds to zero in the
# column rather than to second-guess which small routes are interesting.
ROUTE_SHARE_FLOOR = 0.001

# Clearance kept under the last block on a page, above the page number.
FOOTER_CLEAR = 0.030

# Where the routes table starts, and its line pitch. Both are wanted twice: the
# pass that decides how many pages the routes need has to measure with the same
# numbers the pass that draws them uses, or the last page overruns.
PATHWAY_TOP = 0.868
PATHWAY_PITCH = 0.0128

# How many split channels the isomeric branching block lists. Tungsten has four
# and the page has room; a foil with many would otherwise push the caption off.
SPLITS_SHOWN = 6

# The cover page names the foil. Element symbols only; the benchmark's four
# alloy foils (Inc600, NiCr, SS304, SS316) are already their own names.
# The order the published reports bind a foil's campaigns in, which is also the
# order they are worth reading. 2000exp_5min is the only one covering all 73
# foils, so it leads and is the one a single-experiment report defaults to;
# 1996exp_7hour is last because it is the odd one out, irradiating long enough
# to build up products the 5 minute runs never reach. Anything not named here
# sorts after these, alphabetically, rather than being dropped.
EXPERIMENT_ORDER = ("2000exp_5min", "1996exp_5min", "1996exp_7hour")

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


def half_lives_from(decay=DECAY_LIBRARY, chain=None, cache=CACHE):
    """Half-lives for the ``T1/2`` column, in seconds, keyed by nuclide.

    The one thing this report needs a chain for. Routes and isomeric branching
    come from the solve, filed by step 2, so neither the topology nor the
    branching overlay is read at report time.

    Half-lives live entirely in the decay subsection, so the decay sublibrary
    alone answers this and the library's own chain need not be found.

    ``TransmutationChain`` wants a directory with a manifest and its subsections
    beneath it, and the cache holds the decay subsection's contents rather than a
    chain root, so one symlink and a manifest stand it up. `chain` overrides all
    of it, for a chain that carries its own ``decay/``.

    Returns the mapping and the temporary directory to clean up, or ``({}, None)``
    when there is nothing to read: the column is then left blank, which is what
    it did before and is better than a report that refuses to build over it.
    """
    try:
        import yani

        if chain is not None and (chain / "decay").is_dir():
            print(f"half-lives: {chain}")
            return dict(yani.TransmutationChain(str(chain)).half_lives), None

        borrowed = cache / f"{decay}-transmutation-decay.arrow"
        if not borrowed.is_dir():
            print(f"half-lives: {borrowed} is not there, so the T1/2 column is "
                  f"left blank.\n"
                  f"            Convert the decay sublibrary first, or pass "
                  f"--chain a directory that has its own decay/.")
            return {}, None

        # Held open for as long as the chain is read, and removed after.
        composed = tempfile.TemporaryDirectory(prefix="yani-report-decay-")
        root = pathlib.Path(composed.name)
        (root / "decay").symlink_to(borrowed.resolve(), target_is_directory=True)
        (root / "manifest.json").write_text(json.dumps({
            "format_version": 2,
            "library": decay,
            "subsections": {"decay": {"path": "decay"}},
        }))
        loaded = dict(yani.TransmutationChain(str(root)).half_lives)
        print(f"half-lives: {len(loaded)} from {borrowed.name}")
        return loaded, composed
    except Exception as error:  # noqa: BLE001 - the column is optional
        print(f"half-lives: unreadable ({error}). The T1/2 column is left blank.")
        return {}, None


def edge_weights(result):
    """``(parent, kind, target) -> production per parent atom`` for one result.

    Written by run_transmutation.py from ``get_reaction_rates``. Used to rank
    the isomeric channels by production rather than by rate alone, and empty
    for a decay-only run, which drove no edges.
    """
    out = {}
    for row in result.get("edge_rates") or []:
        parent, kind, target, weight = row
        out[(parent, kind, target)] = weight
    return out


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


def format_percent(value):
    """A sigma as the published tables print it, without rounding a small one to zero.

    Whole percent, because that is what the reports carry and a second decimal
    on a Monte Carlo spread is noise. The exception is a value that rounds to
    zero: ``0%`` and ``-`` would then be the two ways of saying "no uncertainty
    here", one of which means the nuclide is well determined and the other that
    nothing was published about it. ``<1%`` keeps them apart.
    """
    if value is None:
        return "-"
    if 0.0 < value < 0.5:
        return "<1%"
    return f"{value:.0f}%"


def uncertainty_of(result):
    """The nuclear-data sigma on the calculated heat, or None if the run carried none.

    Absent for a run made with ``--no-uncertainty``, and for one made against
    cross sections converted without the covariance. The distinction does not
    matter to a caller: the column is filled when there is a number and left off
    when there is not.
    """
    values = result.get("yani_uncertainty_uW_per_g")
    return None if values is None else np.array(values, dtype=float)


def uncertainty_percent(result):
    """That sigma as a percentage of the value it sits on, per cooling point.

    Relative rather than absolute because that is how the published tables print
    it, and because a percentage is the form that can be compared against the
    measurement's own sigma in the column beside it without dividing anything.
    """
    sigma = uncertainty_of(result)
    if sigma is None:
        return None
    calculated = calculated_of(result)
    return 100.0 * np.divide(sigma, calculated, out=np.zeros_like(sigma),
                             where=calculated > 0)


def coverage_note(result):
    """What the sigma on this page left out, or None if it left nothing out.

    The two uncertainty columns are only honest next to this. A nuclide whose
    evaluation states no MF=33 is perturbed by nothing and comes back exact,
    which is indistinguishable on the page from one that is well determined;
    samples truncated at zero bias the mean upward; and sigmas that had not
    settled are a floor rather than a value. All three are counted by yani
    rather than inferred here, so this cannot claim a coverage the run did not
    have.
    """
    info = result.get("data_uncertainty_info")
    if not info or uncertainty_of(result) is None:
        # No sigma on the page at all, so there are no columns for this to
        # qualify. The caption already says the run carried no covariance;
        # listing every nuclide that had none would say it a second time.
        return None

    parts = []

    # First, because it is the one that decides whether the rest of the sigma
    # is worth reading at all. A library can carry MF=33 for every isotope in
    # the foil and still state nothing about the channel that makes the heat.
    fraction = (info or {}).get("rate_fraction_covered_total")
    if fraction is not None and fraction < COVERAGE_FLOOR:
        parts.append(f"The covariance spans {fraction * 100:.0f}% of the "
                     f"production rate this irradiation drove, so the +/- above "
                     f"is a spread over that part and not over the answer; the "
                     f"channels carrying the rest are perturbed by nothing.")

    missing = info.get("no_covariance_data") or []
    if missing:
        parts.append(f"No MF=33 covariance for {', '.join(sorted(missing))}: "
                     f"reactions on them contribute no sigma, so a product made "
                     f"only from them reads as exact.")

    sampled = info.get("rates_sampled") or 0
    floored = info.get("rates_floored") or 0
    if sampled and floored / sampled > 0.01:
        parts.append(f"{100.0 * floored / sampled:.1f}% of sampled rates went "
                     f"negative and were truncated at zero, biasing the mean up.")

    samples = info.get("samples")
    if samples is not None and not info.get("converged", True):
        parts.append(f"The spread had not settled at {samples} replicas, so these "
                     f"sigmas are a floor.")

    # Deliberately not on the page: `matrices_clipped`, the covariances that were
    # not positive semi-definite and had to be repaired before they could be
    # sampled from. It is real and it is in the JSON and on step 2's stdout, but
    # a repair of a couple of percent is below the precision of every number
    # printed here, and the lines it would cost are the heat curve underneath.
    # What is kept above is the set that changes how a number is read.
    return " ".join(parts) or None


def nuclide_uncertainty_percent(result, nuclide, index):
    """One product's own relative sigma at one cooling point, or None.

    This is the product's share of the calculation moving under the resampled
    cross sections, not its contribution to the total sigma. The two differ
    whenever more than one product carries the heat, so the caption on the page
    says which one it is.
    """
    steps = result.get("by_nuclide_uncertainty_uW_per_g")
    if not steps or index >= len(steps):
        return None
    sigma = steps[index].get(nuclide)
    heat = result["by_nuclide_uW_per_g"][index].get(nuclide)
    if sigma is None or not heat:
        return None
    return 100.0 * sigma / heat


def experiment_order(experiments):
    """The published page order, with anything unrecognised after it by name."""
    return sorted(experiments, key=lambda name: (
        EXPERIMENT_ORDER.index(name) if name in EXPERIMENT_ORDER
        else len(EXPERIMENT_ORDER), name))


def filed_results(results_root):
    """Every ``(foil, experiment)`` a result was filed for, under any library.

    The two roots are the ones ``load_results`` reads, so what this finds is
    exactly what can be reported. Split against the benchmark's own foil names
    rather than on the underscore, because a foil name is not guaranteed to be
    free of one and an experiment name is not guaranteed to be a single field;
    longest foil name first, so ``SS316`` cannot be read as ``SS3``.
    """
    import fns_case

    known = sorted(fns_case.cases(), key=len, reverse=True)
    found = set()
    if not results_root.is_dir():
        return found
    for pattern in ("*/fns_*.json", "*/sweep/fns_*.json"):
        for path in results_root.glob(pattern):
            stem = path.stem[len("fns_"):]
            for case in known:
                if stem.startswith(f"{case}_"):
                    found.add((case, stem[len(case) + 1:]))
                    break
    return found


def discover_cases(results_root):
    """Every foil with a result, so a report can be written for each of them."""
    return sorted({case for case, _experiment in filed_results(results_root)})


def discover_experiments(case, results_root):
    """Every campaign this foil has a result for, in published page order.

    Discovered from the results rather than from the benchmark, because a foil
    measured three times is not thereby a foil that has been run three times,
    and a report naming a page it has no numbers for is worse than one that
    binds what exists and says what is missing.
    """
    return experiment_order({experiment for filed, experiment
                             in filed_results(results_root) if filed == case})


def discover_libraries(case, experiments, results_root):
    """Every library under `results_root` that has a result for this foil.

    Saves naming them: the sweep already filed one folder per library, so the
    report can read off what was run rather than being told again. Ordered by
    name reversed, which puts the later release of a library first and so makes
    it the primary one, whose absolute values the C/E table carries. That is
    what tendl-2025 against tendl-2017 wants, and any pair whose names do not
    sort that way can still be given in the order you want.
    """
    if not results_root.is_dir():
        raise SystemExit(f"no results directory at {results_root}.")
    found = []
    for library in sorted((path.name for path in results_root.iterdir()
                           if path.is_dir()), reverse=True):
        for experiment in experiments:
            for root in (results_root / library / "sweep", results_root / library):
                if (root / f"fns_{case}_{experiment}.json").is_file():
                    found.append(library)
                    break
            else:
                continue
            break
    if not found:
        raise SystemExit(
            f"no {case} results under {results_root}. Run:\n"
            f"  python convert_to_arrow.py --case {case}\n"
            f"  python run_transmutation.py --case {case}")
    return found


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


def nuclide_analysis(result, half_lives, floor=NUCLIDE_EC_FLOOR):
    """The published "nuclide E/C analysis": who carries the heat, and how far off it is.

    One row per product carrying at least `floor` of the calculated heat at some
    one cooling point, read at the point where it carries most.

    The E/C is the total at that point, not the product's own. The measurement
    is a calorimeter reading and does not come apart by nuclide, so there is no
    such thing as one product's E/C, and printing the total against every row
    would be four columns of the same number. What makes a row worth reading is
    the share beside it: an E/C of 0.51 at a point where one product is 98% of
    the calculation is a statement about that product, and the published table
    prints exactly one row for tungsten's 5 minute irradiation for that reason.

    ``%dC_nuc`` is that product's own relative sigma under the resampled
    activation cross sections, read at the same point. It is the column that
    says whether the disagreement beside it is inside what the nuclear data can
    account for, and it is empty for a run made without the covariance.
    """
    values = score(result)
    times = np.array(result["times"], dtype=float)
    rows = []
    for nuclide, share, index in leading_products(result, floor):
        ratio = values["ratio"][index]
        measured = values["measured"][index]
        rows.append({
            "product": nuclide,
            "half_life_s": half_lives.get(nuclide),
            "share_of_calculated": share,
            "at_time": float(times[index]),
            "e_over_c": (1.0 / ratio) if ratio and np.isfinite(ratio) else None,
            "measurement_percent": (100.0 * values["sigma"][index] / measured
                                    if measured else None),
            "calculated_percent": nuclide_uncertainty_percent(result, nuclide, index),
        })
    return rows


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

    For listings that are genuinely preformatted. Anything with columns should
    use :func:`table`, which sets them rather than spacing them.
    """
    font = dict(MONO)
    font.update(style)
    for offset, line in enumerate(lines):
        figure.text(left, top - offset * pitch, line, va="top", **font)
    return top - len(lines) * pitch


def text_height(figure, text, fontsize):
    """Height of a text block as a fraction of the page.

    Measured off the renderer rather than counted in lines times a guessed
    leading. The pathway page divides its height between a table whose length
    the data decides and a note whose wording changes as the code does, so an
    estimate that is 20% out either overruns the page or leaves a third of it
    blank -- both of which this has done.
    """
    artist = figure.text(0.0, 0.0, text, fontsize=fontsize, va="top")
    extent = artist.get_window_extent(figure.canvas.get_renderer())
    artist.remove()
    return extent.height / figure.bbox.height


def table_height(rows, pitch=0.0135, units=False):
    """What :func:`table` will consume from its `top` down to what it returns."""
    return (0.92 * (2 if units else 1) + rows) * pitch + 0.008


def _rule(figure, y, left, right, weight, color=RULE_INK):
    """One horizontal rule across the table's measure."""
    figure.add_artist(plt.Line2D([left, right], [y, y], color=color,
                                 linewidth=weight, transform=figure.transFigure))


def table(figure, columns, rows, top, left=0.08, right=0.94, pitch=0.0135,
          fontsize=7.4, spans=(), rules_after=(), faint_after=()):
    """A set table: booktabs rules, proportional type, figures aligned on their column.

    Columns are `(header, unit, align, weight)`. `weight` divides the measure
    between `left` and `right`, `align` is l/c/r, and `unit` is the second
    header line the published tables carry (uW/g, E/C) or "" for none.

    `spans` are `(label, first, last)`: a heading set across a group of columns
    with a rule under it, for a library that owns both an absolute column and a
    ratio one. `rules_after` are row indices to draw a light rule beneath, for
    separating a summary from the body. `faint_after` are row indices to draw a
    hairline beneath, for a table whose rows come in groups: the pathway table
    gives a product several lines, and without a line closing each group the
    eye loses which product a route belongs to on the way across.

    Numerals are right-aligned on a fixed column edge rather than padded into a
    monospace grid, which is what makes the difference between a typeset table
    and a screenful of terminal output. Returns the height left below.
    """
    total = sum(column[3] for column in columns)
    edges, cursor = [], left
    for _header, _unit, _align, weight in columns:
        width = (right - left) * weight / total
        edges.append((cursor, cursor + width))
        cursor += width

    def place(y, index, text, weight="normal", size=None):
        # A right-aligned column ends one gutter short of its edge, so a figure
        # does not touch the left-aligned column that follows it: without it
        # "1.67m" and "W186(n,2n)W185_m1" set as one word.
        start, end = edges[index]
        align = columns[index][2]
        x = {"l": start, "c": 0.5 * (start + end), "r": end - GUTTER}[align]
        figure.text(x, y, text, ha={"l": "left", "c": "center", "r": "right"}[align],
                    va="top", fontsize=size or fontsize, fontweight=weight)

    # Text is drawn va="top", so a row whose top is at `y` sits in
    # [y - 0.78 * pitch, y]. `above` and `below` put a rule in the gap either
    # side of a row rather than through it.
    above = 0.22 * pitch
    below = 0.86 * pitch

    y = top
    _rule(figure, y + above, left, right, 1.1)

    # A span heads a group of columns, with its own rule showing how far it
    # reaches, in the manner of a cmidrule.
    if spans:
        for label, first, last in spans:
            centre = 0.5 * (edges[first][0] + edges[last][1])
            figure.text(centre, y, label, ha="center", va="top",
                        fontsize=fontsize, fontweight="bold")
        for _label, first, last in spans:
            _rule(figure, y - below, edges[first][0], edges[last][1], 0.5)
        y -= pitch

    for index, (header, _unit, _align, _weight) in enumerate(columns):
        if header:
            place(y, index, header, weight="bold")
    y -= pitch * 0.92

    if any(unit for _h, unit, _a, _w in columns):
        for index, (_header, unit, _align, _weight) in enumerate(columns):
            if unit:
                place(y, index, unit)
        y -= pitch * 0.92

    _rule(figure, y + above, left, right, 0.7)

    for row_index, row in enumerate(rows):
        for index, cell in enumerate(row):
            if cell:
                place(y, index, str(cell))
        y -= pitch
        if row_index in rules_after:
            _rule(figure, y + above, left, right, 0.5)
            y -= pitch * 0.20
        elif row_index in faint_after:
            # No extra space under a hairline: the group rules would otherwise
            # loosen the whole table, and the line alone is enough to close a
            # group. The last group is closed by the bottom rule.
            _rule(figure, y + above, left, right, 0.4, color=FAINT_INK)

    _rule(figure, y + above, left, right, 1.1)
    return y - 0.008


def cover_page(pdf, case, sections, title, subtitle, number, total):
    """The foil, named, and one line per experiment the report covers.

    A foil is measured more than once -- tungsten has two 5 minute irradiations
    at different positions and a 7 hour one -- and the spread between them is a
    result in itself, so the cover carries all of them rather than naming the
    first and leaving the rest to be discovered on page 5.
    """
    figure = page(pdf, title, subtitle, number, total)
    figure.text(0.5, 0.70, element_name(case), ha="center", fontsize=26)
    figure.text(0.5, 0.655, "FNS decay heat validation", ha="center", fontsize=11)

    libraries = ", ".join(library_label(lib) for lib, _ in sections[0][1])
    figure.text(0.5, 0.625, f"libraries: {libraries}", ha="center", fontsize=9,
                color="#555555")

    # The primary library's own uncertainty belongs on the one-glance summary,
    # because a mean deviation is only a verdict next to it: 122% against a
    # calculation good to 6% is a result, and against one good to 60% is not.
    rows = []
    for experiment, results in sections:
        _name, result = results[0]
        values = score(result)
        ratio = np.array(result["ratio"], dtype=float)
        percent = uncertainty_percent(result)
        rows.append([experiment, f"{np.nanmedian(ratio):.3f}",
                     f"{values['mean_percent_diff']:.1f}",
                     "-" if percent is None else f"{np.median(percent):.1f}",
                     f"{values['mean_chi2']:.2f}", f"{len(ratio)}"])
    table(figure, [("Experiment", "", "l", 1.6),
                   ("median C/E", "", "r", 1.0),
                   ("mean % diff", "", "r", 1.0),
                   ("data sigma %", "", "r", 1.1),
                   ("mean chi2", "", "r", 0.9),
                   ("points", "", "r", 0.6)],
          rows, 0.55, left=0.14, right=0.86)
    pdf.savefig(figure)
    plt.close(figure)


# Rows of the volume's ranking that fit on one page, after its heading. Measured
# rather than guessed: the page is A4 at the shared pitch, and a ranking that
# silently truncated at the fold would be worse than one that runs to two pages.
SUMMARY_ROWS = 46


def volume_cover_page(pdf, cases, libraries, title, subtitle, number, total):
    """The volume's own cover: what is bound here, and against what."""
    figure = page(pdf, title, subtitle, number, total)
    figure.text(0.5, 0.66, "FNS decay heat validation", ha="center", fontsize=24)
    figure.text(0.5, 0.605, f"{len(cases)} foils", ha="center", fontsize=13,
                color="#333333")
    figure.text(0.5, 0.565,
                "libraries: " + ", ".join(library_label(l) for l in libraries),
                ha="center", fontsize=9, color="#555555")
    blurb = (
        "One foil says whether a handful of cross sections are right. All of them\n"
        "say whether a library is, which is why they are bound together: the\n"
        "ranking overleaf is the result, and the per-foil pages behind it are\n"
        "where a row that looks wrong is explained.")
    figure.text(0.5, 0.50, blurb, ha="center", va="top", fontsize=9,
                color="#555555", linespacing=1.6)
    pdf.savefig(figure)
    plt.close(figure)


def summary_rows(entries):
    """One ranking row per foil, worst deviation last.

    Ranked on the primary library's mean deviation, and read against the two
    sigmas beside it rather than on its own: a foil 30% out that was measured to
    30% is not a disagreement, and neither is one 30% out whose cross sections
    are known to 30%. Ordered worst last because the point of running every foil
    is to find where the library falls over, and the end of a list is where the
    eye stops.
    """
    rows = []
    for case, sections in entries:
        _experiment, results = sections[0]
        _name, result = results[0]
        values = score(result)
        ratio = np.array(result["ratio"], dtype=float)
        percent = uncertainty_percent(result)
        rows.append((values["mean_percent_diff"], [
            case,
            f"{len(sections)}",
            f"{np.nanmedian(ratio):.3f}",
            f"{values['mean_percent_diff']:.1f}",
            "-" if percent is None else f"{np.median(percent):.1f}",
            f"{result['median_measurement_sigma_percent']:.1f}",
            f"{values['mean_chi2']:.2f}",
        ]))
    return [row for _deviation, row in sorted(rows, key=lambda r: r[0])]


def summary_pages(pdf, entries, title, subtitle, number, total):
    """The ranking, over as many pages as it needs. Returns the next number."""
    columns = [("Foil", "", "l", 1.0),
               ("campaigns", "", "r", 0.9),
               ("median C/E", "", "r", 1.0),
               ("mean % diff", "", "r", 1.0),
               ("data sigma %", "", "r", 1.1),
               ("meas. sigma %", "", "r", 1.2),
               ("mean chi2", "", "r", 0.9)]
    rows = summary_rows(entries)
    chunks = [rows[i:i + SUMMARY_ROWS] for i in range(0, len(rows), SUMMARY_ROWS)] or [[]]
    for index, chunk in enumerate(chunks):
        figure = page(pdf, title, subtitle, number, total)
        heading = "Every foil, ranked by mean deviation"
        if len(chunks) > 1:
            heading += f" ({index + 1} of {len(chunks)})"
        figure.text(0.08, 0.905, heading, fontsize=11)
        bottom = table(figure, columns, chunk, 0.868, left=0.08, right=0.94)
        if index == len(chunks) - 1:
            figure.text(0.08, bottom - 0.020,
                        "Read the deviation against both sigmas beside it. A foil "
                        "out by more than the\nmeasurement can resolve, or than "
                        "its own cross sections are known to, is not a\n"
                        "disagreement the library has to answer for. Each foil's "
                        "own pages follow, in\nthis order.",
                        va="top", fontsize=7.2, color="#555555")
        pdf.savefig(figure)
        plt.close(figure)
        number += 1
    return number


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

    # Laid out as the published tables are: the measurement, then the primary
    # library's own value, then an E/C for every library including that one.
    #
    # A library is named once. The primary's name spans the pair of columns it
    # owns -- its uW/g and its E/C -- rather than sitting above each, which read
    # as two libraries that happened to share a name. Every other library names
    # its single E/C column.
    columns = [("Times", unit[:8], "r", 0.9),
               ("FNS Exp.", "µW/g", "r", 1.5),
               ("", "µW/g", "r", 1.5 if uncertainty_of(primary) is not None else 1.1),
               ("", "E/C", "r", 0.8)]
    columns += [(library_label(name)[:12], "E/C", "r", 0.8) for name, _ in scores[1:]]
    spans = [(library_label(primary_name)[:20], 2, 3)]

    # The calculated column carries its own uncertainty when the run produced
    # one, which is what makes the E/C beside it readable: a ratio of 0.39 is a
    # disagreement if the calculation is good to 6% and barely a statement if it
    # is good to 60%. Without covariance the column is a bare value, as it was.
    calculated_percent = uncertainty_percent(primary)

    rows = []
    for index, moment in enumerate(times):
        measured = primary["measured_uW_per_g"][index]
        sigma = primary["measured_uncertainty"][index]
        percent = 100.0 * sigma / measured if measured else float("nan")
        value = f"{scores[0][1]['calculated'][index]:.2E}"
        if calculated_percent is not None:
            value += f" +/- {format_percent(calculated_percent[index])}"
        row = [f"{moment:.2f}",
               f"{measured:.2E} +/- {percent:.0f}%",
               value]
        for _name, values in scores:
            ratio = values["ratio"][index]
            row.append(f"{1.0 / ratio:.2f}"
                       if ratio and np.isfinite(ratio) else "-")
        rows.append(row)

    blank = [""] * (len(columns) - 2 - len(scores))
    rows.append(["mean % diff. from E", "", *blank]
                + [f"{v['mean_percent_diff']:.0f}" for _, v in scores])
    rows.append(["mean χ²", "", *blank]
                + [f"{v['mean_chi2']:.2f}" for _, v in scores])

    figure.text(0.08, 0.905, f"{element_name(case)}, {primary['experiment']}",
                fontsize=11)
    bottom = table(figure, columns, rows, 0.868, spans=spans,
                   rules_after=(len(times) - 1,))

    # Who carries that heat, so the table above reads as a decay curve rather
    # than as twenty numbers, and so a disagreement has somewhere to land.
    products = leading_products(primary)
    analysis = nuclide_analysis(primary, half_lives)
    figure.text(0.08, bottom - 0.022,
                f"{library_label(primary_name)} nuclide E/C analysis", fontsize=10)
    product_columns = [("Product", "", "l", 1.1),
                       ("T1/2", "", "r", 0.7),
                       ("share of C", "", "r", 0.9),
                       ("at time", unit[:8], "r", 0.8),
                       ("E/C", "", "r", 0.6),
                       ("%ΔE", "", "r", 0.6),
                       ("%ΔCnuc", "", "r", 0.8)]
    product_rows = [[
        row["product"],
        format_half_life(row["half_life_s"]),
        f"{row['share_of_calculated'] * 100:.1f}%",
        f"{row['at_time']:.2f}",
        f"{row['e_over_c']:.2f}" if row["e_over_c"] is not None else "-",
        format_percent(row["measurement_percent"]),
        format_percent(row["calculated_percent"]),
    ] for row in analysis]
    bottom = table(figure, product_columns, product_rows, bottom - 0.040,
                   right=0.72)

    caption = ("E/C is the total at that cooling point, read where the product carries "
               "most of the\ncalculation; a calorimeter reading does not come apart by "
               "nuclide. %ΔE is the\nmeasurement's own sigma there, %ΔCnuc that "
               "product's own spread under its resampled\ncross sections. The +/- on "
               "the µW/g column is the total, taken over whole inventories\nso that a "
               "parent and its daughter move together rather than in quadrature.")
    if uncertainty_of(primary) is None:
        caption = ("E/C is the total at that cooling point, read where the product "
                   "carries most of the\ncalculation; the measurement is a calorimeter "
                   "reading and does not come apart by\nnuclide. %ΔCnuc is empty because "
                   "this run carried no cross-section covariance:\nconvert without "
                   "--no-covariance and run without --no-uncertainty to fill it.")
    # What the sigmas left out, printed under them rather than in a footnote
    # nobody reaches. Absent when the run left nothing out, which is the common
    # case and does not need saying.
    note = coverage_note(primary)
    if note:
        caption += "\n" + "\n".join(textwrap.wrap(note, 86))

    figure.text(0.08, bottom - 0.014, caption, va="top", fontsize=7.2,
                color="#555555")
    bottom -= 0.014 + text_height(figure, caption, 7.2)

    # The primary library's heat curve on the same page as the table it
    # explains, which is how the published pages set it. Sized to whatever the
    # tables left rather than to a fixed box: what is left is not fixed either,
    # and a panel that will not fit is better left to the figure page, which
    # draws it properly, than squeezed over the footer.
    floor = FOOTER_CLEAR + 0.048
    head = bottom - 0.026
    if head - floor >= PANEL_MIN_HEIGHT:
        axes = figure.add_axes([0.115, floor, 0.655, head - floor - 0.016])
        draw_heat(axes, case, primary_name, primary, scores[0][1], products,
                  product_colours(products), fontsize=6.0)

    pdf.savefig(figure)
    plt.close(figure)


def pathway_columns():
    """The routes table's measure, shared by the layout pass and the drawing pass."""
    return [("Product", "", "l", 1.0),
            ("T1/2", "", "r", 0.7),
            ("Pathway", "", "l", 4.2),
            ("path", "", "r", 0.7),
            ("peak share", "", "r", 0.9)]


def pathway_groups(primary, half_lives):
    """One block of rows per product: the line naming it, then its routes.

    Kept apart so a page that cannot hold them all breaks between products
    rather than orphaning a product's routes from the line that names it.

    The routes are read off the result rather than derived here: step 2 asks
    yani, which answers from the chain and the rates the solve actually ran
    with, so the report cannot rank them against a chain that differs from the
    one that produced the numbers.

    Returns the blocks, how many routes were left out of them, and whether the
    result carried routes at all.
    """
    routes = primary.get("production_routes") or {}
    groups, dropped = [], 0

    for nuclide, share, _index in leading_products(primary):
        found = routes.get(nuclide) or []
        carried = len(found)
        # The largest is kept whatever it is, so a product whose every route
        # falls under the floor still names the one that made it rather than
        # going down as a product with no pathway at all.
        above = [r for r in found if r.get("share", 0.0) >= ROUTE_SHARE_FLOOR]
        found = above or found[:1]
        shown = found[:ROUTES_SHOWN]
        dropped += max(0, carried - len(shown))

        if not shown:
            rows = [[nuclide, format_half_life(half_lives.get(nuclide)),
                     "no route within the depth searched", "",
                     f"{share * 100:.1f}%"]]
            groups.append(rows)
            continue

        group = [[nuclide, format_half_life(half_lives.get(nuclide)),
                  shown[0]["route"], f"{shown[0]['share'] * 100:.1f}%",
                  f"{share * 100:.1f}%"]]
        group.extend(["", "", r["route"], f"{r['share'] * 100:.1f}%", ""]
                     for r in shown[1:])
        groups.append(group)

    return groups, dropped, bool(routes)


def isomeric_splits(primary):
    """Every channel that lands in more than one final state, biggest first.

    This is what says whether a disagreement belongs to a cross section or to a
    branching ratio, and it is not in the chain file: the dominant tungsten
    channel carries a placeholder there which the energy-dependent overlay
    replaces at solve time. Step 2 files what the solve used, so this reads it.

    Ranked by production rather than by rate. An edge weight is per atom of its
    parent, so a channel off a 0.12% isotope outranks one off a 28% isotope on
    rate alone, and W180(n,2n) would sit above W186(n,2n) on a tungsten foil.
    Multiplying by the parent's atoms is what makes the order mean anything.
    """
    branching = primary.get("isomeric_branching") or {}
    weights = edge_weights(primary)
    densities = primary.get("initial_atoms_per_barn_cm") or {}
    splits = []
    for parent, kinds in branching.items():
        for kind, split in kinds.items():
            produced = densities.get(parent, 0.0) * sum(
                weight for (a_parent, a_kind, _t), weight in weights.items()
                if (a_parent, a_kind) == (parent, kind))
            splits.append((produced, parent, kind, split))
    return sorted(splits, reverse=True, key=lambda row: row[0])


def isomeric_split_rows(primary):
    """:func:`isomeric_splits` set for the page, capped so the caption still fits."""
    rows = []
    for _produced, parent, kind, split in isomeric_splits(primary)[:SPLITS_SHOWN]:
        # Already ordered largest first by yani.
        shares = "   ".join(f"{target} {fraction * 100:.1f}%"
                            for target, fraction in split)
        rows.append([f"{parent} {kind}", shares])
    return rows


def pathway_preamble(has_routes, dropped):
    """What the routes table means, and what it left out."""
    if not has_routes:
        return ("These results carry no production routes. Step 2 asks yani for "
                "them and files\nthem beside the heat, so re-run it to fill this "
                "page.")

    text = ("Routes are yani's own answer, walked from the nuclides the foil "
            "started with\nover the library's topology and decay data, one "
            "reaction step and then along\nthe decay chain. They are ordered by "
            "\"path\", the share of the product's\nproduction arriving down each "
            "route: the atoms the route starts from, times\nwhat its reaction "
            "drove per atom of its parent over the irradiation, times the\n"
            "branching of every decay it passes through. A route through a 1% "
            "branch\ndelivers 1% of what the reaction made, and a second reaction "
            "step costs another\nfactor of the fluence, which is why one reads "
            "0.0% against a one-step route on\nirradiations this short.\n"
            "\"peak share\" is the product's own largest share of the total decay\n"
            "heat, not a split between its routes.")

    # A cap that is not said out loud reads as "these are all of them". Every
    # product is now carried, so this is only ever about routes.
    if dropped:
        text += (f"\n{dropped} further route{'s' if dropped != 1 else ''} not shown, "
                 f"under {ROUTE_SHARE_FLOOR * 100:.1f}% of their product, "
                 f"or over {ROUTES_SHOWN} per product.")
    return text


def pathway_layout(results, note, half_lives):
    """The whole pathway section, split into pages, before any of it is drawn.

    The footer says "Page n of N", so N has to be known before the first page is
    written, and how many pages the routes need depends on how many rows fit
    under a note whose wording the code decides. Measuring that needs a
    renderer, so a page is made and thrown away here rather than guessed at.

    Every product the C/E table names gets its routes shown, running onto a
    second page where one will not hold them. The truncation this replaces
    dropped whole products for room, and a product with no pathway analysis is
    the one thing the page exists to supply.
    """
    _name, primary = results[0]
    if not primary.get("production_routes"):
        return {"pages": [None], "note": note}

    groups, dropped, has_routes = pathway_groups(primary, half_lives)
    split_rows = isomeric_split_rows(primary)
    preamble = pathway_preamble(has_routes, dropped)

    scratch = plt.figure(figsize=PAGE)
    extras = 0.020 + text_height(scratch, preamble, 7.6)
    if split_rows:
        extras += 0.044 + table_height(len(split_rows))
    plt.close(scratch)

    head = PATHWAY_TOP - table_height(0, PATHWAY_PITCH)
    room_plain = head - FOOTER_CLEAR
    room_last = room_plain - extras

    pages, remaining = [], list(groups)
    while remaining:
        # Everything left fits on a page that also carries the note and the
        # branching block, so this is the last page.
        if sum(len(group) for group in remaining) * PATHWAY_PITCH <= room_last:
            pages.append((remaining, True))
            remaining = []
            break
        taken, used = [], 0
        for group in remaining:
            # The first group goes on whatever its height, so a product with
            # more routes than a page holds still makes progress rather than
            # looping forever on a page it can never fill.
            if taken and (used + len(group)) * PATHWAY_PITCH > room_plain:
                break
            taken.append(group)
            used += len(group)
        # Packing greedily can swallow every product and leave the note a page
        # of its own, which is how 28 tungsten products laid out: one full page
        # and one carrying nothing but the branching block. Hand back as many
        # products as the closing page can hold under the note instead.
        while len(taken) > 1:
            tail = remaining[len(taken) - 1:]
            if sum(len(group) for group in tail) * PATHWAY_PITCH > room_last:
                break
            taken.pop()
        pages.append((taken, False))
        remaining = remaining[len(taken):]
    if not pages or not pages[-1][1]:
        pages.append(([], True))

    return {"pages": pages, "note": None, "preamble": preamble,
            "split_rows": split_rows}


def pathway_page(pdf, case, primary_name, layout, index, title, subtitle,
                 number, total):
    """One page of "how each product that carries heat is made"."""
    figure = page(pdf, title, subtitle, number, total)
    heading = f"{element_name(case)}, {library_label(primary_name)} production pathways"
    if len(layout["pages"]) > 1:
        heading += f" ({index + 1} of {len(layout['pages'])})"
    figure.text(0.08, 0.905, heading, fontsize=11)

    if layout["pages"][index] is None:
        figure.text(0.08, 0.86, layout["note"] or "No routes were available.",
                    va="top", fontsize=8.6, family="monospace")
        pdf.savefig(figure)
        plt.close(figure)
        return

    groups, carries_extras = layout["pages"][index]
    rows, group_ends = [], []
    for group in groups:
        rows.extend(group)
        group_ends.append(len(rows) - 1)

    bottom = PATHWAY_TOP
    if rows:
        bottom = table(figure, pathway_columns(), rows, PATHWAY_TOP,
                       pitch=PATHWAY_PITCH, faint_after=set(group_ends[:-1]))

    if carries_extras:
        if layout["split_rows"]:
            figure.text(0.08, bottom - 0.026,
                        "Flux-weighted isomeric branching, by production", fontsize=10)
            bottom = table(figure, [("Channel", "", "l", 1.0),
                                    ("final states", "", "l", 2.6)],
                           layout["split_rows"], bottom - 0.044, right=0.66)
        figure.text(0.08, bottom - 0.020, layout["preamble"],
                    va="top", fontsize=7.6, color="#555555")

    pdf.savefig(figure)
    plt.close(figure)


def product_colours(products):
    """One colour per product, in the order :func:`leading_products` ranked them.

    tab10 repeats after ten and these pages routinely name more than ten
    products, which put two different nuclides on the same colour. Taking the
    colours from one place is also what keeps a nuclide the same colour in the
    panel under the C/E table and in the panels on the figure page.
    """
    palette = plt.get_cmap("tab20" if len(products) > 10 else "tab10")
    return [palette(index % palette.N) for index in range(len(products))]


def library_label(library):
    """A library's folder name as it should be printed: ``TENDL-2025``.

    These are acronyms, and every published table writes them in capitals. On
    disk they are lowercase because they are directory names, so the two are
    kept apart rather than reconciled: uppercasing here means nothing on disk,
    in the JSON, in the CSV or in ``--libraries`` has to change, and the folder
    a number came from is still the name the page shows.

    ``upper()`` rather than ``title()``, which would give ``Tendl-2025``: the
    letters are an acronym rather than a word, and it is the whole acronym that
    is capitalised. It happens to be right for the release part too, since
    ``endf-b8.1`` wants ``ENDF-B8.1`` and the digits are unaffected.
    """
    return library.upper()


def library_colours(libraries):
    """One colour per library, keyed by name so it is the same on every page.

    Keyed rather than positional because the comparison panel and the figure
    page do not always carry the same libraries in the same order: a library
    with no result for one campaign is dropped from that section only, and a
    colour that shifted along when it did would make the two pages disagree
    about which curve is which.
    """
    palette = plt.get_cmap("tab10")
    return {name: palette(index % palette.N)
            for index, name in enumerate(libraries)}


def draw_band(axes, times, values, sigma, colour):
    """The +/- 1 sigma nuclear-data band around a calculated curve.

    Clipped at the bottom to a floor rather than allowed to reach zero, because
    these panels are logarithmic and a lower edge at or below zero has no place
    on them: matplotlib drops the whole polygon and the band silently vanishes
    at exactly the points where it is widest. The floor is a thousandth of the
    value, which is far enough below the curve to read as "the band runs off the
    bottom here" and cannot be mistaken for a bound.

    A sigma wider than the value itself is not a drawing problem to be hidden:
    it means the cross sections admit an answer near zero, which is a real thing
    for a product made down one poorly known channel.
    """
    lower = np.maximum(values - sigma, values * 1e-3)
    axes.fill_between(times, lower, values + sigma, color=colour,
                      alpha=BAND_ALPHA, linewidth=0, zorder=1)


def heat_limits(values, headroom=PANEL_HEADROOM):
    """Y limits covering the total and the measurement, and deliberately nothing else.

    The panel exists so those two can be read against each other, so those two
    are what set its scale. Ranging on the products instead spends the axis on
    curves that are decades below the total: on tungsten they reach 1e-7 while
    the total never leaves its top three decades, which squeezes the only
    comparison the page is about into the top fifth of the box, and on iron,
    whose total is nearly flat, into a sliver where a 6% disagreement is
    invisible.

    Products that fall below the floor clip, which costs nothing here: they are
    ranked properly on the share panel beside this one, and every one of them is
    in the JSON. The measurement's error bars are included because a lower bar
    can reach under the lowest measured point.

    Returns ``None`` when nothing is positive, which a log axis cannot show
    anyway; the caller then leaves matplotlib's own limits alone.
    """
    edges = np.concatenate([
        np.asarray(values["calculated"], dtype=float),
        np.asarray(values["measured"], dtype=float),
        np.asarray(values["measured"], dtype=float)
        - np.asarray(values["sigma"], dtype=float),
        np.asarray(values["measured"], dtype=float)
        + np.asarray(values["sigma"], dtype=float),
    ])
    edges = edges[np.isfinite(edges) & (edges > 0)]
    if not edges.size:
        return None
    return edges.min() / headroom, edges.max() * headroom


def draw_heat(axes, case, library, result, values, products, colours, fontsize=6.5):
    """The decay-heat curve: each named product, their total, and the measurement."""
    times = np.array(result["times"], dtype=float)
    axes.errorbar(times, values["measured"], yerr=values["sigma"], fmt="^",
                  color="#808080", markersize=3, linewidth=0.7,
                  label="FNS Experiment", zorder=3)
    for index, (nuclide, _share, _at) in enumerate(products):
        axes.plot(times, series_of(result, nuclide), "--", linewidth=0.7,
                  color=colours[index], label=nuclide)
    # The band belongs on the total and on nothing else here. The per-product
    # curves have sigmas too, in the JSON, but drawing thirteen overlapping
    # bands would bury the one comparison the panel is for, and the products'
    # own spreads are already the %dCnuc column on the table page.
    spread = uncertainty_of(result)
    if spread is not None:
        draw_band(axes, times, np.asarray(values["calculated"], dtype=float),
                  spread, "black")
    axes.plot(times, values["calculated"], "-", color="black", linewidth=1.2,
              label="Total", zorder=2)
    axes.set_xscale("log")
    axes.set_yscale("log")
    limits = heat_limits(values)
    if limits is not None:
        axes.set_ylim(*limits)
    axes.set_xlabel(f"Time after irradiation [{result['time_unit']}]", fontsize=fontsize)
    axes.set_ylabel("Heat Output [µW/g]", fontsize=fontsize)
    axes.set_title(f"FNS {result['experiment']} - {case} - {library_label(library)}",
                   fontsize=fontsize + 1.0)
    _finish_panel(axes, fontsize)


def draw_comparison(axes, case, results, colours, fontsize=6.5):
    """Every library's total against the measurement, on one axes.

    The published reports open each foil with this panel, and it is the one
    figure that answers the question the whole document is about: not "is this
    library right" but "do the libraries agree, and is the measurement inside
    them". Four totals drawn separately, one per page, cannot be read against
    each other; drawn together the spread between the evaluations is a distance
    on the page.

    Each total carries its own +/- 1 sigma nuclear-data band, which the
    published figure does not have. Without it the reader has no way to tell a
    gap between two libraries that their own stated uncertainties already cover
    from one that they do not, and those are opposite findings: the first says
    the evaluations are consistent, the second that at least one of them is
    wrong about something it claims to know.

    A library whose evaluations state no MF=33 gets a curve and no band. That is
    the honest rendering, and it is not the same as a narrow band: the caption
    names those libraries so a bare line is not read as a confident one.
    """
    times = np.array(results[0][1]["times"], dtype=float)
    measured = np.array(results[0][1]["measured_uW_per_g"], dtype=float)
    sigma = np.array(results[0][1]["measured_uncertainty"], dtype=float)

    for library, result in results:
        values = calculated_of(result)
        colour = colours[library]
        spread = uncertainty_of(result)
        if spread is not None:
            draw_band(axes, times, values, spread, colour)
        axes.plot(times, values, "-", linewidth=1.1, color=colour,
                  label=library_label(library), zorder=2)

    axes.errorbar(times, measured, yerr=sigma, fmt="^", color="#333333",
                  markersize=3, linewidth=0.7, label="FNS Experiment", zorder=4)

    # Scaled on the totals and the measurement, which is everything drawn here,
    # rather than on `heat_limits`: that one ranges on a single library's
    # calculated array and would cut off whichever other library sits furthest
    # from it, which is the comparison this panel exists to show.
    edges = np.concatenate([measured - sigma, measured + sigma]
                           + [calculated_of(result) for _n, result in results])
    edges = edges[np.isfinite(edges) & (edges > 0)]
    if edges.size:
        axes.set_ylim(edges.min() / PANEL_HEADROOM, edges.max() * PANEL_HEADROOM)

    # Linear in time, and the only panel in this report that is. Every other one
    # is logarithmic because it is showing a decay over three decades, where a
    # log axis is the only one that resolves both ends. This panel is not doing
    # that: it is showing four libraries against one measurement, and it is set
    # the way the published figure sets it so the two can be laid side by side
    # and read as the same picture. A reader checking this report against that
    # one should not have to correct for the axes first.
    #
    # From zero rather than from the first cooling point, again as published.
    axes.set_xlim(left=0.0)
    axes.set_yscale("log")
    axes.set_xlabel(f"Time after irradiation [{results[0][1]['time_unit']}]",
                    fontsize=fontsize)
    axes.set_ylabel("Heat Output [µW/g]", fontsize=fontsize)
    axes.set_title(f"FNS {results[0][1]['experiment']} - {case} - all libraries",
                   fontsize=fontsize + 1.0)
    _finish_panel(axes, fontsize)


def comparison_page(pdf, case, sections, colours, title, subtitle, number, total):
    """One panel per campaign, each carrying every library and the measurement.

    Bound as one page rather than one per campaign, which is how the published
    reports open a foil, and for the same reason the document is one document:
    a foil measured three times is one subject, and the campaigns are worth
    seeing against each other. Tungsten is the case that earns it -- the same
    cross sections read 122%, 65% and 20% out across its three campaigns, and
    three panels on one page say that in a way three pages cannot.
    """
    figure = page(pdf, title, subtitle, number, total)
    rows = len(sections)
    top, floor = 0.90, 0.085

    caption = (
        "Every library's total decay heat against the measurement, one panel per "
        "campaign. The shaded\nband on each curve is that library's own +/- 1 sigma "
        "from its MF=33 cross-section covariance,\nresampled and re-solved per "
        "replica, so two libraries whose bands overlap do not disagree.")
    missing = sorted({library_label(library) for _experiment, results in sections
                      for library, result in results
                      if uncertainty_of(result) is None})
    if missing:
        caption += ("\nDrawn without a band: " + ", ".join(missing) +
                    ". Those runs carried no covariance, which is not the same "
                    "as a\nnarrow one; a bare curve here is an unquantified "
                    "curve, not a confident one.")

    # A band can also be too thin to see for the same reason, which is worse:
    # absence at least looks like absence, while a hairline band looks like a
    # well determined answer. Named here on the same terms, because the panel
    # is where the bands are compared and it is the comparison that misleads.
    thin = sorted({library_label(library) for _experiment, results in sections
                   for library, result in results
                   if ((result.get("data_uncertainty_info") or {})
                       .get("rate_fraction_covered_total") or 1.0) < COVERAGE_FLOOR})
    if thin:
        caption += ("\nBand spans a minority of the production rate: " +
                    ", ".join(thin) + ". Their covariance says nothing about "
                    "the\nchannels carrying most of the heat, so the band is "
                    "narrow because it is nearly empty.")

    # The gridspec bottom is the axes box, and the x label and tick labels of
    # the last panel are drawn below it, so the caption has to clear those as
    # well as itself. XLABEL_CLEAR is that overhang; without it the caption runs
    # through "Time after irradiation".
    height = text_height(figure, caption, 7.2)
    grid = figure.add_gridspec(rows, 1, left=0.10, right=0.82, top=top,
                               bottom=floor + height + XLABEL_CLEAR, hspace=0.42)
    for row, (_experiment, results) in enumerate(sections):
        draw_comparison(figure.add_subplot(grid[row, 0]), case, results, colours)

    figure.text(0.08, floor + height + 0.004, caption, va="top", fontsize=7.2,
                color="#555555")
    pdf.savefig(figure)
    plt.close(figure)


def draw_share(axes, case, library, result, products, colours, fontsize=6.5):
    """The same products as a percentage of the total, which is what a C/E blames."""
    times = np.array(result["times"], dtype=float)
    totals = np.array([sum(step.values()) for step in result["by_nuclide_uW_per_g"]])
    for index, (nuclide, _peak, _at) in enumerate(products):
        percent = np.divide(series_of(result, nuclide) * 100.0, totals,
                            out=np.zeros_like(totals), where=totals > 0)
        axes.plot(times, percent, "--", linewidth=0.7,
                  color=colours[index], label=nuclide)
    axes.set_xscale("log")
    axes.set_ylim(0, 100)
    axes.set_xlabel(f"Time after irradiation [{result['time_unit']}]", fontsize=fontsize)
    axes.set_ylabel("%decay heat contributions", fontsize=fontsize)
    axes.set_title(f"FNS {result['experiment']} - {case} - {library_label(library)}",
                   fontsize=fontsize + 1.0)
    _finish_panel(axes, fontsize)


def _finish_panel(axes, fontsize):
    """The legend outside the axes and the grid, shared by both panel kinds."""
    axes.tick_params(labelsize=fontsize - 0.5)
    axes.legend(fontsize=fontsize - 1.5, loc="upper left", bbox_to_anchor=(1.01, 1.0),
                frameon=False, handlelength=1.5, borderpad=0.2, labelspacing=0.28)
    axes.grid(alpha=0.25, linewidth=0.4)


def figure_page(pdf, case, results, title, subtitle, number, total):
    """Heat curves left, % contributions right, one row per library."""
    figure = page(pdf, title, subtitle, number, total)
    rows = len(results)
    top, floor = 0.90, 0.085

    # One library gets the two panels stacked at full width, which is the shape
    # the published pages use and the only one that fills an A4 page. Side by
    # side, two panels 2.6 inches wide cannot: laying a single row across the
    # full height drew each panel nearly four times taller than wide and
    # squashed a curve spanning three decades into a vertical smear, and sizing
    # it by aspect instead left two thirds of the page blank.
    #
    # More than one library keeps a row per library, panels side by side, sized
    # by aspect so the rows shrink to fit rather than stretching. Two panels
    # abreast fix the width, so few libraries leave height over; the block is
    # centred in what is left rather than hung from the top, which put all of
    # the slack in one hole under the last row. `bottom` still clears the
    # footer: with enough libraries they compress instead of running over it.
    if rows == 1:
        grid = figure.add_gridspec(2, 1, left=0.10, right=0.82, top=top,
                                   bottom=0.16, hspace=0.30)
        cells = [(grid[0, 0], grid[1, 0])]
    else:
        left, right = 0.09, 0.90
        hspace = 0.42
        panel_w = (right - left) / (2 + PANEL_WSPACE)
        stack = rows + hspace * (rows - 1)
        # The aspect that would fill the page exactly, held between the drawn
        # shape and the tallest one worth setting.
        fitted = (top - floor) / (stack * panel_w * PAGE[0] / PAGE[1])
        aspect = min(PANEL_ASPECT_MAX, max(PANEL_ASPECT, fitted))
        panel_h = aspect * panel_w * PAGE[0] / PAGE[1]
        span = panel_h * stack
        slack = max(0.0, (top - floor) - span)
        grid = figure.add_gridspec(rows, 2, left=left, right=right,
                                   top=top - slack / 2,
                                   bottom=max(floor, top - slack / 2 - span),
                                   hspace=hspace, wspace=PANEL_WSPACE)
        cells = [(grid[row, 0], grid[row, 1]) for row in range(rows)]

    for row, (library, result) in enumerate(results):
        heat_cell, share_cell = cells[row]
        products = leading_products(result)
        colours = product_colours(products)
        draw_heat(figure.add_subplot(heat_cell), case, library, result,
                  score(result), products, colours)
        draw_share(figure.add_subplot(share_cell), case, library, result,
                   products, colours)

    pdf.savefig(figure)
    plt.close(figure)


def pathway_data(primary, half_lives):
    """The pathway analysis as data: every product, every route, every share.

    The page has to fit on a page, so it caps the routes per product and drops
    the ones under :data:`ROUTE_SHARE_FLOOR`. This does neither. A sidecar that
    reproduced the page's elisions would be no more use than the page, and the
    route somebody wants to check is exactly the one small enough to have been
    cut.
    """
    routes = primary.get("production_routes") or {}
    return [{
        "product": nuclide,
        "half_life_s": half_lives.get(nuclide),
        "peak_share": share,
        "routes": [{"route": r["route"], "path": r["share"],
                    "atoms_per_barn_cm": r["production"]}
                   for r in routes.get(nuclide) or []],
    } for nuclide, share, _index in leading_products(primary)]


def write_data(out, case, sections, half_lives):
    """Every table the PDF draws, beside it, in a form something else can read.

    A PDF is where this report is read and a poor place to get a number back
    out of. The JSON carries the lot, uncapped and unrounded; the CSV carries
    just the C/E table, which is the one that goes into a spreadsheet.

    Written unconditionally rather than behind a flag. The whole report takes
    under a second to regenerate, so a sidecar nobody asked for costs nothing,
    and a flag people have to know about is the thing that makes an output hard
    to use.
    """
    document = {
        "case": case,
        "element": element_name(case),
        "experiments": [],
    }
    for experiment, results in sections:
        _primary_name, primary = results[0]
        document["experiments"].append({
            "experiment": experiment,
            "time_unit": primary["time_unit"],
            "times": primary["times"],
            "measured_uW_per_g": primary["measured_uW_per_g"],
            "measured_uncertainty_uW_per_g": primary["measured_uncertainty"],
            "libraries": [{
                "library": library,
                "calculated_uW_per_g": calculated_of(result).tolist(),
                "calculated_uncertainty_uW_per_g": (
                    None if uncertainty_of(result) is None
                    else uncertainty_of(result).tolist()),
                "e_over_c": [(1.0 / r) if r and np.isfinite(r) else None
                             for r in score(result)["ratio"]],
                "mean_percent_diff": score(result)["mean_percent_diff"],
                "mean_chi2": score(result)["mean_chi2"],
                "data_uncertainty_info": result.get("data_uncertainty_info"),
            } for library, result in results],
            "nuclide_analysis": nuclide_analysis(primary, half_lives),
            "pathways": pathway_data(primary, half_lives),
            "isomeric_branching": [
                {"parent": parent, "reaction": kind, "final_states": split}
                for _produced, parent, kind, split in isomeric_splits(primary)],
        })

    data_path = out.with_suffix(".json")
    data_path.write_text(json.dumps(document, indent=2) + "\n")

    csv_path = out.with_name(out.stem + "_ce.csv")
    libraries = [library for library, _ in sections[0][1]]
    header = ["experiment", "time", "time_unit", "measured_uW_per_g",
              "measured_uncertainty_uW_per_g"]
    for library in libraries:
        header += [f"{library}_uW_per_g", f"{library}_uncertainty_uW_per_g",
                   f"{library}_E_over_C"]
    with csv_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for block in document["experiments"]:
            for index, moment in enumerate(block["times"]):
                row = [block["experiment"], moment, block["time_unit"],
                       block["measured_uW_per_g"][index],
                       block["measured_uncertainty_uW_per_g"][index]]
                for entry in block["libraries"]:
                    sigma = entry["calculated_uncertainty_uW_per_g"]
                    row += [entry["calculated_uW_per_g"][index],
                            "" if sigma is None else sigma[index],
                            entry["e_over_c"][index]]
                writer.writerow(row)

    return [data_path, csv_path]


def prepare(case, experiments, libraries, results_root, half_lives):
    """Everything a foil's pages need, and how many pages that will be.

    Separated from the drawing because the footer says "Page n of N" and N is
    not known until the routes have been laid out. A volume needs every foil's
    count before it can draw the first one, so the counting happens here and the
    drawing in :func:`draw_case`.
    """
    sections = [(experiment, load_results(case, experiment, libraries, results_root))
                for experiment in experiments]

    # The routes come with each result, because they are a statement about the
    # rates that result was solved with as much as about the topology.
    note = None
    if not any(result.get("production_routes")
               for _experiment, results in sections for _library, result in results):
        note = ("These results carry no production routes. Re-run step 2 to "
                "fill this page.")

    layouts = [pathway_layout(results, note, half_lives)
               for _experiment, results in sections]
    pages = 1 + sum(2 + len(layout["pages"]) for layout in layouts)
    return sections, layouts, pages


def draw_case(pdf, case, sections, layouts, colours, half_lives,
              title, subtitle, number, total):
    """One foil's pages, into an already-open PDF. Returns the next page number.

    The comparison panel first, then per campaign the C/E table, the pathways
    and the figures. Takes the page number rather than starting at one, so the
    same pages sit in a single-foil report or partway through a volume without
    knowing which they are in.
    """
    comparison_page(pdf, case, sections, colours, title, subtitle, number, total)
    number += 1
    for (_experiment, results), layout in zip(sections, layouts):
        table_page(pdf, case, results, half_lives, title, subtitle, number, total)
        number += 1
        for index in range(len(layout["pages"])):
            pathway_page(pdf, case, results[0][0], layout, index,
                         title, subtitle, number, total)
            number += 1
        figure_page(pdf, case, results, title, subtitle, number, total)
        number += 1
    return number


def report_scores(case, sections):
    """The per-campaign figures of merit, printed as the report is written."""
    for experiment, results in sections:
        for library, result in results:
            values = score(result)
            print(f"{case:>8} {experiment:>16} {library:>12}: "
                  f"mean % diff {values['mean_percent_diff']:6.1f}   "
                  f"mean chi^2 {values['mean_chi2']:8.2f}   "
                  f"{len(leading_products(result))} products above "
                  f"{SHARE_FLOOR * 100:g}% peak share")


def build(case, experiments, libraries, results_root, chain, decay, out,
          title, subtitle):
    """Write one foil's report and return where it went.

    One cover, then the foil's own pages. A foil measured more than once belongs
    in one document: the spread between its experiments is a result, and the
    published report binds them the same way.
    """
    half_lives, composed = half_lives_from(decay, chain)
    sections, layouts, pages = prepare(case, experiments, libraries,
                                       results_root, half_lives)
    filed = sections[0][1][0][1].get("production_routes") or {}
    print(f"pathways: {sum(len(v) for v in filed.values())} routes into "
          f"{len(filed)} products")
    total = 1 + pages
    colours = library_colours(libraries)

    out.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(out) as pdf:
        cover_page(pdf, case, sections, title, subtitle, 1, total)
        draw_case(pdf, case, sections, layouts, colours, half_lives,
                  title, subtitle, 2, total)
        pdf.infodict()["Title"] = f"{title}: {element_name(case)}"
        pdf.infodict()["Subject"] = subtitle

    if composed is not None:
        composed.cleanup()

    written = write_data(out, case, sections, half_lives)
    report_scores(case, sections)
    return [out, *written]


def build_volume(cases, libraries, results_root, chain, decay, out,
                 title, subtitle, experiments=None):
    """Every foil in one document: cover, ranking, then each foil's pages.

    73 separate PDFs are 73 documents nobody opens. What running every foil buys
    is the comparison between them, and that only exists as a page if they are
    bound together. The per-foil pages are unchanged and follow the ranking, so
    a row that looks wrong is a page turn away from the pathways that explain it.

    A foil with no results is skipped and named rather than failing the volume:
    on a sweep of 73 there will usually be a few, and losing the other 70 to one
    of them would be absurd.
    """
    half_lives, composed = half_lives_from(decay, chain)

    prepared, skipped = [], []
    for case in cases:
        found = experiments or discover_experiments(case, results_root)
        if not found:
            skipped.append(case)
            continue
        try:
            sections, layouts, pages = prepare(case, found, libraries,
                                               results_root, half_lives)
        except SystemExit as error:
            print(f"{case}: skipped ({error})")
            skipped.append(case)
            continue
        prepared.append((case, sections, layouts, pages))

    if not prepared:
        raise SystemExit(f"no foil under {results_root} has results to bind.")

    entries = [(case, sections) for case, sections, _l, _p in prepared]
    ranking = max(1, -(-len(entries) // SUMMARY_ROWS))
    total = 1 + ranking + sum(pages for _c, _s, _l, pages in prepared)
    colours = library_colours(libraries)

    out.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(out) as pdf:
        volume_cover_page(pdf, [c for c, *_ in prepared], libraries,
                          title, subtitle, 1, total)
        number = summary_pages(pdf, entries, title, subtitle, 2, total)
        # Bound in the ranking's order, so the volume reads the way the summary
        # points: the foils a reader was sent to looking for are together.
        order = {row[0]: index for index, row in enumerate(summary_rows(entries))}
        for case, sections, layouts, _pages in sorted(
                prepared, key=lambda item: order.get(item[0], 0)):
            number = draw_case(pdf, case, sections, layouts, colours, half_lives,
                               title, subtitle, number, total)
        pdf.infodict()["Title"] = f"{title}: {len(prepared)} foils"
        pdf.infodict()["Subject"] = subtitle

    if composed is not None:
        composed.cleanup()
    for case, sections, _l, _p in prepared:
        report_scores(case, sections)
    if skipped:
        print(f"skipped {len(skipped)}: {' '.join(skipped)}")
    print(f"\n{len(prepared)} foils, {total} pages")
    return [out]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--case", nargs="+", default=None,
                        help="foils to report, as named in the benchmark (W, Fe, "
                             "SS316), one document each. Default: every foil that "
                             "has a result under --results")
    parser.add_argument("--experiments", nargs="+", default=None,
                        help="which experiments to report, in page order. Default: "
                             "every campaign the foil has a result for, which is "
                             "what puts one foil's whole measurement history in one "
                             "document: tungsten has 2000exp_5min, 1996exp_5min and "
                             "1996exp_7hour, and the spread between them is a result")
    parser.add_argument("--libraries", nargs="+", default=None,
                        help="folders under --results, one per library. The first "
                             "is the primary one, whose absolute values the table "
                             "carries. Default: every folder that has a result for "
                             "this foil, newest-looking library first")
    parser.add_argument("--results", type=pathlib.Path, default=HERE / "results",
                        help="root the sweeps were filed under (default: results/)")
    parser.add_argument("--volume", action="store_true",
                        help="bind every selected foil into one document: a "
                             "cover, a ranking of all of them, then each foil's "
                             "own pages in that order. What running every foil "
                             "buys is the comparison between them, and that is "
                             "only a page if they are bound together")
    parser.add_argument("--chain", type=pathlib.Path, default=None,
                        help="a chain directory carrying its own decay/, to take "
                             "the T1/2 column from instead of --decay. Only "
                             "half-lives are read from it: the routes and the "
                             "isomeric branching come with the results")
    parser.add_argument("--decay", default=DECAY_LIBRARY,
                        help="decay sublibrary supplying half-lives when --chain has "
                             f"no decay/ of its own (default: {DECAY_LIBRARY})")
    parser.add_argument("--output", type=pathlib.Path, default=None,
                        help="PDF to write, which only makes sense for a single "
                             "foil (default: results/report_<case>.pdf, or "
                             "results/report_<case>_<experiment>.pdf when one "
                             "experiment was named)")
    parser.add_argument("--title", default="YANI decay heat validation",
                        help="running header, top line")
    parser.add_argument("--subtitle", default="FNS benchmark, JAEA / IAEA CoNDERC",
                        help="running header, second line")
    args = parser.parse_args()

    cases = args.case or discover_cases(args.results)
    if not cases:
        raise SystemExit(
            f"no results under {args.results}. Run:\n"
            f"  python convert_to_arrow.py --case Fe\n"
            f"  python run_transmutation.py --case Fe")
    if args.output is not None and len(cases) > 1 and not args.volume:
        raise SystemExit(
            f"--output names one file but {len(cases)} foils were selected "
            f"({', '.join(cases)}). Name one foil with --case, bind them with "
            f"--volume, or drop --output and they will be written as "
            f"results/report_<case>.pdf.")
    if args.case is None:
        print(f"foils: {', '.join(cases)} (found under {args.results})")

    if args.volume:
        libraries = args.libraries or discover_libraries(
            cases[0], discover_experiments(cases[0], args.results), args.results)
        out = args.output or (args.results / "report_volume.pdf")
        for path in build_volume(cases, libraries, args.results, args.chain,
                                 args.decay, out, args.title, args.subtitle,
                                 args.experiments):
            print(f"\n  pdf: {path}")
        return

    for case in cases:
        # One document per foil, carrying every campaign that foil was measured
        # in. A foil measured more than once is one subject, and the spread
        # between its campaigns is a result about the data rather than about any
        # one measurement: iron reads 6% high against 2000exp_5min and 7% low
        # against 1996exp_5min, which is only visible when they are bound
        # together.
        experiments = args.experiments or discover_experiments(case, args.results)
        if not experiments:
            print(f"{case}: no results under {args.results}, skipped")
            continue
        if args.experiments is None:
            print(f"{case}: {', '.join(experiments)}")

        libraries = args.libraries or discover_libraries(case, experiments,
                                                         args.results)
        if not args.libraries:
            print(f"libraries: {', '.join(libraries)}, primary {libraries[0]} "
                  f"(found under {args.results})")

        # Named for the foil, because that is what the document is about. The
        # experiment only reaches the filename when one was asked for by name,
        # where it is the whole of what distinguishes the file.
        stem = (f"{case}_{experiments[0]}" if args.experiments
                and len(args.experiments) == 1 else case)
        out = args.output or (args.results / f"report_{stem}.pdf")
        written = build(case, experiments, libraries, args.results,
                        args.chain, args.decay, out, args.title,
                        args.subtitle)
        print()
        for path in written:
            print(f"{path.suffix.lstrip('.'):>6}: {path}")
        print()


if __name__ == "__main__":
    main()
