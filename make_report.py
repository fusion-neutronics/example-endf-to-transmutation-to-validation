#!/usr/bin/env python3
"""Step 5: bind the sweep results into a decay-heat validation report PDF.

Reads the per-foil JSON that ``run_transmutation.py`` and ``sweep_fns.py``
already write, one sweep per library, and lays them out in the page order the
UKAEA decay-heat validation reports use:

    make_report.py                 # one document per foil, all of its campaigns
    make_report.py --case W        # just this one

    cover page               the foil, named, and one line per campaign
    C/E table                one row per cooling point, one E/C column per
                             library, the nuclide E/C analysis under it saying
                             which product carries the disagreement, and that
                             library's heat curve under that
    production pathways      the reaction and decay steps that make each
                             product, over as many pages as they need
    figure page              heat curves and % contributions, one row per library

The last three repeat per campaign, because one foil is one document. A foil
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

Every column the published reports carry is now filled. Three were structural
omissions rather than oversights, and are recorded here because the reason each
one was empty says what it is:

* An uncertainty on the calculated value (the reports print "+/- 6%"), and with
  it ``%dC_nuc`` in the nuclide analysis. Both are cross-section covariance
  carried through the inventory. yani-core 0.11.0 reads the MF=33 covariance
  step 1 now writes beside each nuclide, resamples the activation cross sections
  from it, folds each draw against the foil's own spectrum and re-solves, so the
  spread over that ensemble is the number. Step 2 writes it into its JSON and
  both columns are read from there.

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

* Decay feeds. ``TransmutationChain`` gained ``decays`` in yani 0.8, so a route
  can be followed past its neutron reaction and ``W186(n,2n)W185_m1(IT)W185`` is
  derivable rather than guessed at.
* "Path %", the fraction of a product's inventory arriving down each route, and
  with it the flux-weighted isomeric branching. The solve computed the rate of
  every edge to build its burnup matrix and then discarded it, so the weights
  could not be recovered from the result. yani-core 0.9.0 hands them back
  (``TransmutationResults.get_reaction_rates``, fusion-neutronics/core#505),
  run_transmutation.py writes them into its JSON, and both are now on the page.

  The branching is the number that says which evaluated quantity a disagreement
  belongs to, and it cannot be read off the chain: the chain file's own
  branching for the dominant tungsten channel is a placeholder,

      W186 (n,2n) -> W185      branching 1.000000
      W186 (n,2n) -> W185_m1   branching 0.000000

  which the energy-dependent overlay in ``branching/`` replaces at solve time.
  It is not zero in any meaningful sense: with the overlay configured W185m
  carries 98% of the decay heat at the first cooling point, and without it 10%.
  The split the solve actually used is now printed rather than inferred.

That is why this script needs yani 0.11.1 or newer. It still reads a result
filed by an older run: one with no per-edge rates falls back to the unweighted
route order and says so rather than printing zeros, and one with no sigma leaves
the two uncertainty columns empty rather than at zero. What it will not do is
invent either from what is there.

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

# How many decay steps a route may follow after its neutron reaction, and the
# branching below which a decay branch is not worth a line. Three steps covers
# what these irradiations reach: the longest route the published pages carry is
# W182(n,p)Ta182n(IT)Ta182m(IT)Ta182, which is three.
DECAY_STEPS = 3
DECAY_FLOOR = 0.01

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
#
# Only applied when the result carries per-edge rates. Without them there is no
# share to compare against and every route stands, which is the same fallback
# the ordering takes.
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


def route_source(chain, seeds, depth, case, decay_floor=DECAY_FLOOR,
                 decay_steps=DECAY_STEPS):
    """Production routes per product, as whole strings like the reports print them.

    ``W186(n,2n)W185_m1(IT)W185``: one neutron reaction off an isotope of the
    foil, then the decay steps that carry the product on. Both halves come from
    the chain, :attr:`reactions` for the first and :attr:`decays` for the rest.

    Two scopes are applied to the reaction half, and both matter.

    The chain covers every nuclide the library knows, so asking it what makes
    W187 answers with Os190(n,a) and Ir192(n,npa) as readily as W186(n,gamma).
    Those are real edges and irrelevant here: nothing in a tungsten foil is
    osmium. So the walk starts at the foil's own isotopes.

    Walking from there to saturation does not fix it, because successive capture
    does reach osmium and iridium from tungsten eventually. It takes a fluence
    the FNS foils never see: these irradiations are 5 minutes to 7 hours, where
    the inventory is first and second order in the reaction rate and nothing
    else registers. ``depth`` is the cap that keeps the page to that, and it
    counts reaction steps only. Decay does not spend it, because a decay is not
    a fluence-dependent step: Hf183 becomes Ta183 on its own schedule whether
    the beam is on or not.

    Decay is capped separately by ``decay_steps``, and branches below
    ``decay_floor`` are dropped. Cycles are refused outright, which is not
    hypothetical: Ta186 beta-decays to W186, an isotope of the foil, so an
    unguarded walk goes round for as long as it is allowed to.
    """
    if chain is None:
        return {}, ("No chain was given, so the routes are unavailable. Pass "
                    "--chain the\ndirectory convert_to_arrow.py wrote for this "
                    "foil to fill this page.")
    forward = chain.reactions
    # New in yani 0.8: without it a walk stops at the first product, and the
    # decay steps that carry it on are exactly what a pathway analysis is made
    # of. See the module docstring.
    decays = chain.decays

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

    routes = {}

    def record(product, text, origin, steps, edges):
        routes.setdefault(product, {})[text] = (origin, steps, edges)

    def follow_decays(node, text, origin, steps, visited, budget, edges):
        """Extend one route along the decay chain, depth first."""
        for mode, daughter, branching in decays.get(node, []):
            if not daughter or daughter in visited or branching < decay_floor:
                continue
            extended = f"{text}({mode}){daughter}"
            record(daughter, extended, origin, steps + 1, edges)
            if budget > 1:
                follow_decays(daughter, extended, origin, steps + 1,
                              visited | {daughter}, budget - 1, edges)

    # Breadth first over the reaction steps, carrying the route text so a
    # second reaction extends the first rather than replacing it.
    frontier = [(seed, seed, seed, {seed}, ()) for seed in covered]
    for _ in range(depth):
        following = []
        for node, text, origin, visited, edges in frontier:
            for reaction, target, _ratio in forward.get(node, []):
                if not target or target in visited:
                    continue
                extended = f"{text}{reaction}{target}"
                steps = extended.count("(")
                # Every reaction edge the route crosses, in order. All of them
                # rather than the first: a route is only as likely as its least
                # likely step, and pricing it on one edge ranked a triple
                # capture above the direct channel. The decay steps that follow
                # carry the same edges -- what reached the product came down
                # these reactions and then decayed.
                walked = edges + ((node, reaction, target),)
                record(target, extended, origin, steps, walked)
                follow_decays(target, extended, origin, steps,
                              visited | {target}, decay_steps, walked)
                following.append((target, extended, origin,
                                  visited | {target}, walked))
        frontier = following
    return routes, None


def edge_weights(result):
    """``(parent, kind, target) -> production per parent atom`` for one result.

    Written by run_transmutation.py from ``get_reaction_rates``, which arrived
    with yani-core 0.9.0. Absent from results filed before that, and absent for
    a decay-only run, in which case the routes fall back to being unweighted.
    """
    out = {}
    for row in result.get("edge_rates") or []:
        parent, kind, target, weight = row
        out[(parent, kind, target)] = weight
    return out


def route_shares(product_routes, weights, densities):
    """Fraction of a product's production arriving down each of its routes.

    The reports call this "Path %". An edge weight is the reactions that edge
    drove per atom of its parent over the whole irradiation, so it is
    dimensionless and small -- of order 1e-4 on these foils. A route's share is
    the atoms it starts from times that weight once per reaction step it
    crosses.

    Multiplying along the route is what keeps a multi-step route in its place.
    Pricing one edge and calling the rest free made
    ``W184(n,gamma)W185(n,gamma)W186(n,2n)W185_m1`` outrank the direct
    ``W186(n,2n)W185_m1``, which is wrong by eight orders of magnitude: each
    further reaction costs another factor of the fluence, and these
    irradiations are 5 minutes to 7 hours. Only TENDL's topology showed it --
    endf-b8.1 does not carry the intermediate edges, so every route it offered
    was one step and the two prices agreed.

    It is an ordering, not an inventory: the exact k-step term carries a
    combinatorial factor this leaves out. That does not reach the ranking,
    which the fluence decides.

    Decay steps do not enter. A route's decay tail moves what the reactions
    made, it does not make more of it, so two routes differing only past the
    last reaction share its production rather than splitting it.

    Returns ``{}`` when nothing can be weighted, which is the honest answer for
    a result with no rates in it: the caller keeps the unweighted order and the
    page says the column is unavailable rather than printing zeros.
    """
    scored = {}
    for text, (origin, _steps, edges) in product_routes.items():
        weight = densities.get(origin, 0.0)
        for edge in edges:
            weight *= weights.get(edge, 0.0)
        if weight > 0.0:
            scored[text] = weight
    total = sum(scored.values())
    if total <= 0.0:
        return {}
    return {text: weight / total for text, weight in scored.items()}


def isomeric_split(weights, parent, kind):
    """The flux-weighted branching of one channel over its final states.

    The split arrives as two edges of one channel, so the sum is the channel's
    production and the ratio is what the reports call the isomeric branching.
    It cannot be read off the chain file: its own branching for the dominant
    tungsten channel is a placeholder, ``W186 (n,2n) -> W185`` at 1.0 and
    ``-> W185_m1`` at 0.0, which the energy-dependent overlay replaces at solve
    time. The difference between the two is a factor of 50 in the answer.

    Returns ``{target: fraction}``, or ``{}`` when the channel drove nothing.
    """
    targets = {target: weight
               for (a_parent, a_kind, target), weight in weights.items()
               if a_parent == parent and a_kind == kind and target is not None}
    total = sum(targets.values())
    if total <= 0.0:
        return {}
    return {target: weight / total for target, weight in sorted(targets.items())}


def route_order(item):
    """Sort key for the routes into one product: fewest steps first, then by name.

    Used when the result carries no per-edge rates, which is every result filed
    before yani-core 0.9.0. It carries no claim beyond short routes before long
    ones, and the page says so.

    Weighting it by anything guessed is worse than not weighting it. An earlier
    version ordered by the natural abundance of the isotope a route starts from,
    which put ``W184(n,gamma)W185_m1`` above ``W186(n,2n)W185_m1`` because W184
    is 30.6% against W186's 28.4%. That is backwards: in a 14 MeV spectrum
    capture is negligible against (n,2n), and W185m is 98% of tungsten's decay
    heat, so the row a reader looks at first was the row the heuristic got
    wrong. Abundance was standing in for the rate; now the rate is available,
    :func:`route_shares` uses it and this is the fallback.
    """
    text, (_origin, steps, _edges) = item
    return steps, text


def shared_route_order(shares):
    """Sort key ordering routes by share, biggest first, then as :func:`route_order`."""
    def key(item):
        text, (_origin, steps, _edges) = item
        return -shares.get(text, 0.0), steps, text
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

    Absent for a run made with ``--no-uncertainty``, for one made against cross
    sections converted without the covariance, and for every result filed before
    yani-core 0.11.0, which could not produce it. The distinction does not matter
    to a caller: the column is filled when there is a number and left off when
    there is not, which is what it did for every library before this existed.
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


def discover_chain(libraries, data_root):
    """The chain convert_to_arrow.py wrote beside one of these libraries' cross sections.

    The primary library's own is preferred, and that is not a nicety: its
    ``branching/`` is the energy-dependent overlay, and without it the dominant
    tungsten channel takes the chain file's placeholder and W185m goes from 98%
    of the decay heat to 10%.

    Falling back to another library on the report is still worth doing. A result
    can outlive the conversion it came from, and the pathway page mostly asks
    the chain for topology and half-lives, which the libraries agree on far
    better than they agree on the branching. Returns which library it came from
    so the caller can say, because a chain quietly borrowed from the library the
    page is not about is exactly the kind of thing that should be said out loud.
    """
    for library in libraries:
        chain = data_root / library / "chain"
        if chain.is_dir():
            return chain, library
    return None, None


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

    libraries = ", ".join(lib for lib, _ in sections[0][1])
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
    columns += [(name[:12], "E/C", "r", 0.8) for name, _ in scores[1:]]
    spans = [(primary_name[:20], 2, 3)]

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
                f"{primary_name} nuclide E/C analysis", fontsize=10)
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


def pathway_groups(primary, routes, half_lives):
    """One block of rows per product: the line naming it, then its routes.

    Kept apart so a page that cannot hold them all breaks between products
    rather than orphaning a product's routes from the line that names it.
    Returns the blocks, how many routes were left out of them, and whether the
    result carried the rates needed to rank routes at all.
    """
    weights = edge_weights(primary)
    densities = primary.get("initial_atoms_per_barn_cm") or {}
    groups, dropped, weighted_any = [], 0, False

    for nuclide, share, _index in leading_products(primary):
        product_routes = routes.get(nuclide, {})
        shares = route_shares(product_routes, weights, densities)
        weighted_any = weighted_any or bool(shares)
        order = shared_route_order(shares) if shares else route_order
        found = sorted(product_routes.items(), key=order)
        carried = len(found)
        if shares:
            # The largest is kept whatever it is, so a product whose every route
            # falls under the floor still names the one that made it rather than
            # going down as a product with no pathway at all.
            above = [item for item in found
                     if shares.get(item[0], 0.0) >= ROUTE_SHARE_FLOOR]
            found = above or found[:1]
        shown = [text for text, _meta in found[:ROUTES_SHOWN]]
        dropped += max(0, carried - len(shown))
        if not shown:
            shown = ["no route within the depth searched"]

        def path_column(text):
            """The route's share, blank when this product could not be weighted."""
            return f"{shares[text] * 100:.1f}%" if text in shares else ""

        group = [[nuclide, format_half_life(half_lives.get(nuclide)),
                  shown[0], path_column(shown[0]), f"{share * 100:.1f}%"]]
        group.extend(["", "", text, path_column(text), ""] for text in shown[1:])
        groups.append(group)

    return groups, dropped, weighted_any


def isomeric_splits(primary):
    """Every channel that lands in more than one final state, biggest first.

    This is what says whether a disagreement belongs to a cross section or to a
    branching ratio, and it is not in the chain file: the dominant tungsten
    channel carries a placeholder there which the energy-dependent overlay
    replaces at solve time.

    Ranked by production rather than by rate. An edge weight is per atom of its
    parent, so a channel off a 0.12% isotope outranks one off a 28% isotope on
    rate alone, and W180(n,2n) would sit above W186(n,2n) on a tungsten foil.
    Multiplying by the parent's atoms is what makes the order mean anything.
    Every channel is returned rather than only the largest, because the largest
    by production is not always the one that carries the heat.
    """
    weights = edge_weights(primary)
    densities = primary.get("initial_atoms_per_barn_cm") or {}
    splits = []
    for parent, kind in {(parent, kind) for parent, kind, _target in weights}:
        split = isomeric_split(weights, parent, kind)
        if len(split) < 2:
            continue
        produced = densities.get(parent, 0.0) * sum(
            weight for (a_parent, a_kind, _t), weight in weights.items()
            if (a_parent, a_kind) == (parent, kind))
        splits.append((produced, parent, kind, split))
    return sorted(splits, reverse=True, key=lambda row: row[0])


def isomeric_split_rows(primary):
    """:func:`isomeric_splits` set for the page, capped so the caption still fits."""
    rows = []
    for _produced, parent, kind, split in isomeric_splits(primary)[:SPLITS_SHOWN]:
        shares = "   ".join(f"{target} {fraction * 100:.1f}%"
                            for target, fraction in sorted(split.items(),
                                                           key=lambda kv: -kv[1]))
        rows.append([f"{parent} {kind}", shares])
    return rows


def pathway_preamble(depth, weighted_any, dropped):
    """What the routes table means, and what it left out."""
    if weighted_any:
        provenance = (
            "ordered by \"path\", the share of the product's production arriving down\n"
            "each route. That is the atoms the route starts from times, once per\n"
            "reaction step it crosses, the reactions that step drove per atom of its\n"
            "own parent -- from the per-edge rates the solve computed and now hands\n"
            "back. Each further reaction step therefore costs another factor of the\n"
            "fluence, which is why a two-step route reads 0.0% against a one-step one\n"
            "on irradiations this short. Decay steps do not enter: a decay moves what\n"
            "the reactions made rather than making more of it. A blank \"path\" is a\n"
            "route whose reactions drove nothing the solve could price.\n")
    else:
        provenance = (
            "ordered fewest steps first and then by name, which is not a ranking of\n"
            "any kind. \"path\" needs the per-edge reaction rates, which this result\n"
            "was filed without: re-run it against yani 0.9.0 or newer to fill the\n"
            "column.\n")
    text = ("Routes are the library's own topology and decay data, walked "
            f"{depth} reaction step{'s' if depth != 1 else ''}\n"
            "from the foil's own isotopes and then along the decay chain. They are\n"
            + provenance +
            "\"peak share\" is the product's own largest share of the total decay\n"
            "heat, not a split between its routes.")

    # A cap that is not said out loud reads as "these are all of them". Every
    # product is now carried, so this is only ever about routes.
    if dropped:
        limit = (f"under {ROUTE_SHARE_FLOOR * 100:.1f}% of their product, "
                 f"or over {ROUTES_SHOWN} per product" if weighted_any
                 else f"at most {ROUTES_SHOWN} per product")
        text += (f"\n{dropped} further route{'s' if dropped != 1 else ''} not shown, "
                 f"{limit}.")
    return text


def pathway_layout(results, routes, note, half_lives, depth):
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
    if not routes:
        return {"pages": [None], "note": note}

    groups, dropped, weighted_any = pathway_groups(primary, routes, half_lives)
    split_rows = isomeric_split_rows(primary)
    preamble = pathway_preamble(depth, weighted_any, dropped)

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
    heading = f"{element_name(case)}, {primary_name} production pathways"
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
    axes.plot(times, values["calculated"], "-", color="black", linewidth=1.2,
              label="Total", zorder=2)
    axes.set_xscale("log")
    axes.set_yscale("log")
    limits = heat_limits(values)
    if limits is not None:
        axes.set_ylim(*limits)
    axes.set_xlabel(f"Time after irradiation [{result['time_unit']}]", fontsize=fontsize)
    axes.set_ylabel("Heat Output [µW/g]", fontsize=fontsize)
    axes.set_title(f"FNS {result['experiment']} - {case} - {library}",
                   fontsize=fontsize + 1.0)
    _finish_panel(axes, fontsize)


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
    axes.set_title(f"FNS {result['experiment']} - {case} - {library}",
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


def pathway_data(primary, routes, half_lives):
    """The pathway analysis as data: every product, every route, every share.

    The page has to fit on a page, so it caps the routes per product and drops
    the ones under :data:`ROUTE_SHARE_FLOOR`. This does neither. A sidecar that
    reproduced the page's elisions would be no more use than the page, and the
    route somebody wants to check is exactly the one small enough to have been
    cut.
    """
    weights = edge_weights(primary)
    densities = primary.get("initial_atoms_per_barn_cm") or {}
    products = []
    for nuclide, share, _index in leading_products(primary):
        product_routes = routes.get(nuclide, {})
        shares = route_shares(product_routes, weights, densities)
        order = shared_route_order(shares) if shares else route_order
        products.append({
            "product": nuclide,
            "half_life_s": half_lives.get(nuclide),
            "peak_share": share,
            "routes": [{"route": text, "path": shares.get(text)}
                       for text, _meta in sorted(product_routes.items(), key=order)],
        })
    return products


def write_data(out, case, sections, half_lives, routes, depth):
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
        "pathway_depth": depth,
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
            "pathways": pathway_data(primary, routes, half_lives),
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


def build(case, experiments, libraries, results_root, chain, decay, depth, out,
          title, subtitle):
    """Write the report and return where it went.

    One cover, then three pages per experiment. A foil measured more than once
    belongs in one document: the spread between its experiments is a result, and
    the published report binds them the same way.
    """
    sections = [(experiment, load_results(case, experiment, libraries, results_root))
                for experiment in experiments]
    loaded, composed = read_chain(chain, decay)
    half_lives = dict(loaded.half_lives) if loaded else {}

    # The routes are the chain's own topology walked from the foil's isotopes,
    # so they do not depend on which irradiation was run. The weights on them
    # do: each experiment sits at its own position in its own spectrum, and the
    # per-edge rates come from that experiment's own result.
    routes, note = route_source(loaded, foil_isotopes(case, experiments[0]),
                                depth, case)
    if note:
        print(f"pathways: {note.splitlines()[0]}")
    else:
        print(f"pathways: {sum(len(v) for v in routes.values())} routes into "
              f"{len(routes)} products")

    # Laid out before any page is drawn, because the footer says "of N" and the
    # routes decide how many pages they need.
    layouts = [pathway_layout(results, routes, note, half_lives, depth)
               for _experiment, results in sections]
    total = 1 + sum(2 + len(layout["pages"]) for layout in layouts)

    out.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(out) as pdf:
        cover_page(pdf, case, sections, title, subtitle, 1, total)
        number = 2
        for (_experiment, results), layout in zip(sections, layouts):
            table_page(pdf, case, results, half_lives, title, subtitle,
                       number, total)
            number += 1
            for index in range(len(layout["pages"])):
                pathway_page(pdf, case, results[0][0], layout, index,
                             title, subtitle, number, total)
                number += 1
            figure_page(pdf, case, results, title, subtitle, number, total)
            number += 1
        pdf.infodict()["Title"] = f"{title}: {element_name(case)}"
        pdf.infodict()["Subject"] = subtitle

    if composed is not None:
        composed.cleanup()

    written = write_data(out, case, sections, half_lives, routes, depth)

    for experiment, results in sections:
        for library, result in results:
            values = score(result)
            print(f"{experiment:>16} {library:>12}: "
                  f"mean % diff {values['mean_percent_diff']:6.1f}   "
                  f"mean chi^2 {values['mean_chi2']:8.2f}   "
                  f"{len(leading_products(result))} products above "
                  f"{SHARE_FLOOR * 100:g}% peak share")
    return [out, *written]


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
    parser.add_argument("--chain", type=pathlib.Path, default=None,
                        help="chain directory convert_to_arrow.py wrote, for "
                             "half-lives and the pathway page. Its branching/ is "
                             "what puts isomers on the page, so prefer the "
                             "library's own over a generic chain. Default: "
                             "data/<primary library>/chain if it is there, and "
                             "without it those columns are left out")
    parser.add_argument("--data", type=pathlib.Path, default=HERE / "data",
                        help="root convert_to_arrow.py filed its conversions "
                             "under, searched for the default --chain "
                             "(default: data/)")
    parser.add_argument("--decay", default=DECAY_LIBRARY,
                        help="decay sublibrary supplying half-lives when --chain has "
                             f"no decay/ of its own (default: {DECAY_LIBRARY})")
    parser.add_argument("--pathway-depth", type=int, default=1,
                        help="how many reaction steps from the foil's own isotopes "
                             "the pathway page walks, not counting decay steps, "
                             "which are followed regardless. 1 matches the published "
                             "pages, whose routes are all one reaction and then "
                             "decay; 2 adds second-order routes like "
                             "W186(n,3n)W184(n,gamma)W185m, which need a fluence "
                             "these 5 minute to 7 hour irradiations never reach "
                             "(default: 1)")
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
    if args.output is not None and len(cases) > 1:
        raise SystemExit(
            f"--output names one file but {len(cases)} foils were selected "
            f"({', '.join(cases)}). Name one foil with --case, or drop --output "
            f"and they will be written as results/report_<case>.pdf.")
    if args.case is None:
        print(f"foils: {', '.join(cases)} (found under {args.results})")

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

        chain = args.chain
        if chain is None:
            chain, from_library = discover_chain(libraries, args.data)
            if chain is not None and from_library != libraries[0]:
                print(f"chain: taken from {from_library}, which is not the primary "
                      f"library. Its isomeric branching is {from_library}'s, not "
                      f"{libraries[0]}'s")

        # Named for the foil, because that is what the document is about. The
        # experiment only reaches the filename when one was asked for by name,
        # where it is the whole of what distinguishes the file.
        stem = (f"{case}_{experiments[0]}" if args.experiments
                and len(args.experiments) == 1 else case)
        out = args.output or (args.results / f"report_{stem}.pdf")
        written = build(case, experiments, libraries, args.results,
                        chain, args.decay, args.pathway_depth, out, args.title,
                        args.subtitle)
        print()
        for path in written:
            print(f"{path.suffix.lstrip('.'):>6}: {path}")
        print()


if __name__ == "__main__":
    main()
