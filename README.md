# ENDF to transmutation to validation

An end-to-end example: raw evaluated nuclear data goes in one end, and a
comparison against a measured experiment comes out the other.

The whole chain is three steps.

1. **ENDF in.** Evaluations for the natural isotopes of whatever the foil is
   made of. TENDL-2025 by default, any library you point it at otherwise.
2. **Arrow out.** NJOY reconstructs the resonances and Doppler broadens to
   294 K, and `nuclear_data_to_arrow` writes the result as Arrow, which is the
   form yani reads. This happens on your machine; nothing pre-processed is
   downloaded.
3. **Validation.** yani irradiates a 1 g foil with the measured neutron
   spectrum from the JAEA FNS decay-heat experiment, follows the activation
   products through the cooling schedule, and the specific decay heat is
   plotted against the measurement, with an uncertainty band on each: the
   experiment's own, and the one the cross-section covariance puts on the
   calculation.

It defaults to iron, and `--case` takes any of the 73 FNS foils.

There is no neutron transport anywhere in this example. The FNS experiment
published the spectrum the foil actually saw, so that spectrum is the input and
yani collapses it against the cross sections to get one-group reaction rates.
That is what yani is for: transmutation and activation given a spectrum.

## Install

### Clone the repo
```bash
git clone https://github.com/fusion-neutronics/example-endf-to-transmutation-to-validation.git
cd example-endf-to-transmutation-to-validation
```

### Create an isolated Python environment
```bash
python3 -m venv .venv
source .venv/bin/activate
```

### Install the dependencies needed by the example
```bash
python -m pip install -r requirements.txt
```

That is the whole install: everything comes from PyPI, nothing from git, and
nothing needs a checkout beside this one. `yani` both reads the nuclear data
and makes it: the ENDF conversion, the reaction topology and the isomeric
branching are all functions on the wheel.

The wheels are abi3 from Python 3.10 up, one per platform rather than one per
Python version, published for Linux x86_64, Windows x64 and macOS arm64. There
is no Rust toolchain to install and nothing is compiled here.

NJOY is the one external tool left, and no library in any language avoids it:
an ENDF evaluation describes the resonance region with parameters rather than
pointwise cross sections, so something has to reconstruct and Doppler broaden
it. `requirements.txt` pulls the `njoy2016` wheel, which puts an `njoy` in the
venv. Use a build of your own with `--njoy`, which is also how to point at the
IAEA-NDS fork that FENDL needs; both builds succeed, so the wrong choice is
quiet.

## Run

```bash
python convert_to_arrow.py
python run_transmutation.py
```

That is iron, the default. Both scripts take `--case` for any of the 73 FNS
foils, and both need the same one, since step 1 converts exactly the isotopes
step 2 will ask for:

```bash
python run_transmutation.py --list      # every foil and its experiments

python convert_to_arrow.py  --case W
python run_transmutation.py --case W

python convert_to_arrow.py  --case SS316    # an alloy: 6 elements, 24 isotopes
python run_transmutation.py --case SS316
```

Compounds and alloys need no special handling. The foil's composition is read
from the benchmark's own input deck, so `--case Cs` picks up the caesium
carbonate the experiment actually used (Cs 81.58%, O 14.73%, C 3.69%) and
converts the isotopes of all three elements.

### Which experiment

Each foil was measured in up to three campaigns, chosen with `--experiment`:

| experiment | foils | irradiation | measured over | median sigma |
|---|---|---|---|---|
| `2000exp_5min` (default) | 73 | 5 min | 0.6 to 55 minutes | 6.8% |
| `1996exp_5min` | 30 | 5 min | 1.1 to 57 minutes | 5.0% |
| `1996exp_7hour` | 29 | 7.6 hours | 0.6 to 400 days | 6.1% |

`2000exp_5min` is the default because it is the only one that covers every
foil, and it is what `openmc_activator`'s own comparison notebook uses. The
1996 5-minute set has tighter uncertainties on the foils it does cover. The
7-hour set is the interesting one physically: it irradiates long enough and
follows the foil far enough out to test long-lived products that a 5-minute
irradiation never builds up.

### Pointing it at your own ENDF files

`convert_to_arrow.py` needs an ENDF evaluation for each natural isotope of each
element in the foil. Give it a directory holding them and it will find them:

```bash
python convert_to_arrow.py --endf-dir /path/to/endf --library endf-b8.1
```

The directory is searched recursively, and every file in it is identified by
the nuclide named in its own header rather than by its name. Naming conventions
do not matter: `n-Fe056.tendl`, `n-026_Fe_056.endf` and `whatever.dat` are all
recognised, non-ENDF files in the same directory are ignored, and metastable
evaluations are excluded by their LIS0 flag rather than by spotting an `m` in
the filename. Scanning a 2850-file library takes under a second. If an isotope
is missing, or two evaluations of the same one are present, it says which and
stops rather than guessing.

`--library` is the name stamped into the Arrow output, so set it to whatever
you are actually converting. `$ENDF_DIR` works in place of `--endf-dir`.

With no `--endf-dir` it fetches the library `--library` names, for you:

```bash
python convert_to_arrow.py                          # TENDL-2025, the default
python convert_to_arrow.py --library tendl-2017     # or TENDL-2017
```

The two archives are laid out differently inside, TENDL-2025 flat at the root
and TENDL-2017 nested under `neutron_file/<El>/<El><A>/lib/endf/`, which costs
you nothing: evaluations are matched on their filename wherever they sit.

The last form transfers around 3 GB whatever you asked for, because TENDL
publishes its neutron sublibrary as a single archive; only the wanted members
are written to disk. Conversion runs about twenty seconds per isotope, and each
run reconverts the isotopes for the selected foil.

`run_transmutation.py` writes `results/<source>/fns_<case>_<experiment>.png` and
a JSON of the same numbers, including the per-nuclide breakdown of the heat and
the nuclear-data sigma on both.

When step 1 uses `--endf-dir`, the `<source>` name defaults to that directory's
basename. For example, this pair matches:

```bash
python convert_to_arrow.py --case W --endf-dir libs --library tendl-2025
python run_transmutation.py --case W --source libs
```

or, equivalently:

```bash
python run_transmutation.py --case W --cross-sections data/libs/neutron
```

When `--cross-sections` ends with `/neutron`, `run_transmutation.py` defaults
`--chain` to the sibling `/chain-<case>` directory.

The foil is in that name because the chain is scoped to one foil's isotopes
while the cross sections are not. `data/tendl-2025/neutron/` accumulates every
nuclide ever converted into it, so converting tungsten after iron adds
`W186.arrow` beside `Fe56.arrow` and takes nothing away; a chain written to a
shared path would instead replace the previous foil's, leaving that foil
unrunnable and the report walking its pathways through a topology built for
somebody else's parents.

### One folder per library

Everything a run produces is filed under the nuclear data it came from:
`data/tendl-2017/`, `results/tendl-2025/`, and so on. Evaluations of your own
are filed under the directory they came from, and `--source` overrides the name
when the directory is not what you want it called:

```bash
python sweep_fns.py --library tendl-2017 --endf-dir data/tendl-2017-endf
python sweep_fns.py --library tendl-2025 --source tendl-2025 \
    --endf-dir /path/to/a/directory/called/something/else
```

This is not tidiness. TENDL names its evaluations identically in every release,
so `n-Fe056.tendl` from 2017 and from 2025 are the same filename. Unpacked into
one directory the second library silently reuses the first one's files, and a
run you believe is TENDL-2017 reports TENDL-2025 numbers.

### Comparing two libraries

```bash
python compare_libraries.py tendl-2017 tendl-2025
```

Pairs the two sweeps foil by foil on `|C/E - 1|`. A foil only one library
managed is excluded rather than counted as a loss for the other, and a
difference smaller than that foil's measurement sigma is called a tie, because
most of these foils agree within it and a difference below sigma is not a
result.

The tie threshold is the measurement's sigma alone, which is now the narrower of
the two available. `sweep_fns.py` carries a `data sigma` column per foil, the
median nuclear-data uncertainty on that foil's calculated heat, and on several
foils it is the larger number: a 6% difference between two libraries on a foil
whose cross sections are each known to 6% is not a result either. Folding it in
would need a view on how correlated two evaluations of the same physics are,
which is a real question and not one this repo answers, so the threshold is left
where the published comparisons put it and the second sigma is printed beside it
to be read.

### A report PDF

```bash
python make_report.py            # one document per foil, all of its campaigns
python make_report.py --case W   # just this one
```

Binds the per-foil JSON into `results/report_<case>.pdf`: the foil named; then
every library's total against the measurement, one panel per campaign; then the
C/E table with an E/C column per library and, under it, the nuclide E/C analysis
and that library's heat curve; then the production pathways, over as many pages
as they need; then the heat curves with their percentage contributions. Nothing
is recomputed: step 2 already wrote the full per-nuclide breakdown and the sigma
on it, so a report is cheap to regenerate and cannot disagree with the run it
came from.

**One document per foil, covering every campaign it was measured in.** The last
three pages repeat per campaign. A foil measured more than once is still one
subject, and the spread between its campaigns is a result about the data rather
than about any one measurement, so binding them separately hides it. Tungsten's
three campaigns run 122%, 65% and 20% out, all against the same cross sections;
iron reads 6% high on `2000exp_5min` and 7% low on `1996exp_5min`.

With no arguments it writes one for every foil that has a result. Which foils,
which campaigns, which libraries and which chain are all read off what is on
disk. `--case` picks foils, `--experiments` picks campaigns, and naming exactly
one campaign puts it back in the filename as `report_<case>_<experiment>.pdf`.

The calculation's own uncertainty is on the page beside its value, `+/- 6%` and
`%ΔCnuc` both. A result filed before those existed is still readable: the
two uncertainty columns are left empty rather than filled with a zero, which
would be a different claim.

`--libraries` is also how to say which library is the primary one, whose
absolute values the table carries:

```bash
python make_report.py --case W --libraries tendl-2025 tendl-2017
```

Beside the PDF go the same tables in a form something else can read:
`report_<case>.json`, which carries every table for every campaign with none of
the page's caps applied, and `report_<case>_ce.csv`, which carries the C/E
tables alone, one row per campaign per cooling point. A PDF is where this report
is read and a poor place to get a number back out of.

The table's two figures of merit are the conventional ones, the mean of
`|C/E - 1|` and a mean chi-squared against the measurement's own sigma, which is
what says whether a deviation is larger than the experiment can resolve. On
W/2000exp_5min against TENDL-2025 they come out at 122.0% and 105.14. Neither
carries the calculation's own uncertainty, because neither does by definition;
that number is in the column beside the value it belongs to instead.

The report needs no chain. Routes and isomeric branching come with the results,
and the only thing left that a chain could answer is the `T1/2` column, which is
half-lives and so lives entirely in the decay sublibrary the run already used.
That is read straight from the cache, so there is no library chain to find and
none to go stale: a report whose chain has moved on is no longer a report whose
pathway page quietly changes.

`--chain` remains for a chain that carries its own `decay/`, and `--decay` names
the sublibrary otherwise. Neither is needed for the common case.

Routes are whole strings, `W186(n,2n)W185m(IT)W185`, one neutron reaction off an
isotope of the foil and then the decay steps that carry the product on. **yani
works them out, not this repo.** Step 2 asks
`TransmutationResults.get_production_routes` and files the answer beside the
heat; step 5 binds it. That matters because a route is a statement about the
rates a solve ran with as much as about the topology, and deriving it later
means loading the chain a second time and hoping it is the one that produced the
numbers.

They are walked from the foil's own isotopes rather than across the whole chain:
asking the chain what makes W187 answers with `Os190(n,a)` and `Ir192(n,npa)` as
readily as `W186(n,gamma)`, and nothing in a tungsten foil is osmium. One
reaction step, which is what the published pages carry, and then up to three
decay steps.

"Path %" is the share of a product's production arriving down each route: the
atoms the route starts from, times what its reaction drove per atom of its
parent over the irradiation, times the branching of every decay it passes
through. A route through a 1% decay branch delivers 1% of what the reaction
made. Each further reaction step carries the step duration too, so a two-step
route is in the same units as a one-step one and comes out smaller by roughly a
factor of the fluence, which on a 5 minute irradiation at 1e10 n/cm2/s is parts
in a billion.

The **flux-weighted isomeric branching** comes from the same place,
`get_isomeric_branching`, and is the number that says whether a disagreement
belongs to a cross section or to a branching ratio. It is not in the chain file:
the dominant tungsten channel carries a placeholder `0.0` there, which the
energy-dependent overlay replaces at solve time with the 53.5/46.5 the spectrum
actually gives.

Routes carrying under 0.1% of their product are left off the page and counted
in the note under it. A route that made a thousandth of a product says nothing
about where a disagreement lives, and six lines of `0.0%` under every product
made down one channel read as if the column were not a ranking at all. The JSON
beside the PDF carries every route, weighted, with nothing dropped.

The **nuclide E/C analysis** under the C/E table names the products that could
account for a disagreement: each one's largest share of the calculation, when
it reaches it, and the E/C there. The E/C is the total at that cooling point,
not the product's own, because the measurement is a calorimeter reading and
does not come apart by nuclide. What makes the row worth reading is the share
beside it: an E/C of 0.51 at a point where one product is 98% of the
calculation is a statement about that product.

### The uncertainty on the calculated value

The last column of that table, `%ΔCnuc`, and the `+/- 6%` the published reports
print beside every calculated µW/g, are both cross-section covariance carried
through the inventory.

Step 1 writes the MF=33 covariance as a `covariance.arrow` beside each nuclide.
That is on by default, because it is nearly free: on W186 it costs no measurable
NJOY time and 78 kB against a 2.0 MB nuclide. Step 2 then resamples the
activation cross sections from it, folds each draw against the foil's own
spectrum and re-solves the schedule, and the spread over that ensemble is the
number. The solver is untouched and only its input moves, so the mean inventory
is what it always was: `--no-uncertainty` reproduces the old bare value exactly,
and the whole ensemble costs about a second.

On W/2000exp_5min against TENDL-2025 the calculated heat comes out good to 5.5%
in the median, 4.3% to 6.1% over the cooling points, and `%ΔCnuc` on W185m comes
out at 6%. Both are printed beside the value they qualify, so a deviation is
never read without them.

Two things about that number are worth knowing before it is used.

It is **the activation cross sections**, in practice and not by construction.
Half-lives, decay branching ratios, fission yields and the isomeric-branching
overlay are held at their evaluated values. The flux is the one that could move
and does not: yani perturbs the spectrum as a second source and asks for
both by default (`DataUncertainty.available_sources()` returns
`['cross_sections', 'flux_spectrum']`), but the FNS benchmark publishes its
measured spectrum as 709 group fluxes with no sigma on them, so there is nothing
to draw from and the source contributes exactly zero. yani counts that rather
than passing over it: `spectra_without_flux_sigma` is 1 on every run here.

The distinction matters because those are different statements. A flux held
fixed by choice would be an approximation this example made; a flux held fixed
for want of a published sigma is a gap in the benchmark, and it will close on
its own the day the spectrum arrives with uncertainties. Step 2 prints what it
did not propagate rather than leaving that to be assumed.

A **sigma of zero is not a claim of certainty**. An evaluation that states no
MF=33 is perturbed by nothing and its products come back exact, which on the
page is indistinguishable from a product that is well determined. yani counts
those nuclides, step 2 names them, and the report page prints a note under the
table rather than letting the column read as coverage it does not have. The
tables also print `<1%` rather than `0%` for a small non-zero spread, so that
`0%` is never the answer to two different questions.

The total is taken over **whole inventories** rather than assembled from
per-nuclide sigmas. Decay heat is a function of an inventory, and a parent and
its daughter move together under one resampled cross section, so adding
per-nuclide sigmas in quadrature would assume an independence the resampling
exists precisely to avoid assuming. `%ΔCnuc` is the narrower thing: one
product's own spread, which is what says whether the row beside it is inside
what that product's cross sections can account for.

## What comes out

For iron on the default experiment, yani tracks the measurement to 6.0% mean
deviation, starting on top of it and drifting high by the end of the hour, with
a median C/E of 1.064 against a measurement whose own median sigma is 5.4%. Run
the same iron against `1996exp_5min` and it comes out low instead, median C/E
0.927 for 7.3%. The two campaigns bracket the calculation, which is a fair
reminder that a single C/E number is a statement about one measurement and not
only about the data.

The breakdown says why the shape is what it is. Mn56 from Fe56(n,p), half-life
2.58 hours, is 78% of the heat at the first point and 99.9% at the last. The
rest of the early heat is Mn57 and Fe53, both of which are gone within minutes,
so the first few points test three cross sections and the tail tests one.

### What the sigma is a statement about

The uncertainty on the calculated value is the library's own account of itself,
not a property of the calculation, and iron is the case that makes that
unmissable. The same foil, the same spectrum, the same solver, on
`2000exp_5min`:

| library | median C/E | mean deviation | data sigma | isotopes with MF=33 | rate covered |
|---|---|---|---|---|---|
| tendl-2025 | 1.064 | 6.0% | **33%** | 4 of 4 | 99.9% |
| jeff-4.0 | 1.060 | 5.5% | **1.2%** | 4 of 4 | 99.5% |
| endf-b8.1 | 1.057 | 4.9% | **1.2%** | 2 of 4 | 95.6% |
| jendl-5.0 | 1.094 | 8.8% | **4.9%** | 2 of 4 | 89.2% |

All four track the measurement to within about 9%, and the sigmas beside that
span a factor of 27. Nothing in the transport, the chain or the solver accounts
for it: it is what each evaluation states about its own Fe56(n,p), which carries
most of this foil's heat. TENDL's covariances come from varying TALYS parameters
and are wide; ENDF/B-VIII.1's Fe56 is one of the most measured cross sections
there is and its covariance is correspondingly tight.

Read the sigma against the deviation and they say different things: a 6%
deviation against a 33% sigma is a foil that library cannot resolve, and the
same 6% against a genuine 1.2% would be a real disagreement.

**The last column is the one that decides whether the sigma beside it means
anything**, and it is not the isotope count. Counting isotopes with MF=33 asks
whether an evaluation says something; weighting by production asks whether it
says it about the reactions this irradiation actually drove.

On iron the count is the pessimistic one. ENDF/B-VIII.1 covers 2 of 4 isotopes,
which reads as half a foil, and the two it covers are Fe54 and Fe56; Fe56 alone
is 91.75% of natural iron, so 95.6% of the production is covered and the 1.2%
beside it is a real number. JENDL-5 is the same shape. Reading the isotope count
alone would throw away two sigmas that are worth having.

Tungsten is where it fails the other way, and completely:

| library | median C/E | mean deviation | data sigma | isotopes with MF=33 | rate covered |
|---|---|---|---|---|---|
| tendl-2025 | 2.355 | 122.0% | 5.5% | 5 of 5 | 99.3% |
| jeff-4.0 | 1.805 | 78.9% | <0.1% | 5 of 5 | **5.6%** |
| endf-b8.1 | 1.805 | 78.9% | <0.1% | 5 of 5 | **5.6%** |
| jendl-5.0 | 1.669 | 71.7% | none | 0 of 5 | none |

JEFF-4.0 and ENDF/B-VIII.1 state covariance for every natural tungsten isotope,
so "5 of 5" is true and reads as complete coverage. What they state it for is
`(n,3n)` and `(n,gamma)`. The `W186(n,2n)W185m` that makes 98% of this foil's
decay heat has no MF=33 at all, so the ensemble perturbs 5.6% of the production
and reports a spread of under 0.1%: the most confident number in the table and
the least earned. yani files the per-channel coverage as
`rate_fraction_covered`, step 2 folds it against its own reaction rates and the
foil's own densities, and the report page says so under any table whose
covariance spans less than 90% of the production.

Both weights are load-bearing. Weighting by rate alone, without the density of
the parent each rate belongs to, counts a channel on a trace isotope the same as
one on the bulk and puts ENDF/B-VIII.1's iron at 32% rather than 96%, which is
the opposite conclusion about whether its sigma is worth reading.

JENDL-5 is the honest end of the same problem. Its tungsten carries no MF=33
anywhere, so no ensemble runs, and the columns are left empty rather than
filled with a zero. An empty column is a question; a `0.0%` is an answer, and it
would be the wrong one.

On tungsten against TENDL-2025 the calculated heat is good to 5.5% in the median
while sitting 122% from the measurement, which is a disagreement no account of
the cross-section uncertainty absorbs, and the pathway page says where it lives.

Swapping the cross sections for another library, with everything else held
fixed, turns this into a library comparison:

```bash
python convert_to_arrow.py --endf-dir /path/to/endfb-8.1 --library endf-b8.1 \
    --output data/b81
python run_transmutation.py --cross-sections data/b81 --output results/b81
```

### Four libraries at once

The same foil against as many libraries as you have, each converted from its own
ENDF on your own machine:

```bash
for lib in tendl-2025 jeff-4.0 endf-b8.1 jendl-5.0; do
  python convert_to_arrow.py  --case W --endf-dir /path/to/$lib --source $lib
  python run_transmutation.py --case W --source $lib
done
python make_report.py --case W --libraries tendl-2025 jeff-4.0 endf-b8.1 jendl-5.0
```

`--libraries` fixes the column order and which library is primary; without it the
report picks its own. Tungsten on `2000exp_5min`, against a measurement whose own
median sigma is 13.1%:

| library | median C/E | mean deviation | mean chi2 | data sigma |
|---|---|---|---|---|
| tendl-2025 | 2.355 | 122.0% | 105.14 | 5.5% |
| jeff-4.0 | 1.805 | 78.9% | 42.81 | <0.1% |
| endf-b8.1 | 1.805 | 78.9% | 42.81 | <0.1% |
| jendl-5.0 | 1.669 | 71.7% | 37.05 | none |

Every one of them is further from the measurement than any of them is from the
others, which is the useful shape of that table: this is a disagreement about
tungsten, not about which library to pick.

One row is worth reading twice. JEFF-4.0 and ENDF/B-VIII.1 come out **identical**,
to the digit, and that is a fact about the evaluations rather than a coincidence:
JEFF-4.0 adopted ENDF/B-VIII.1's tungsten wholesale. The two files on disk agree
byte for byte over MF=2, MF=3, MF=8 and MF=33. What differs is MF=1, which is
descriptive text, and one endpoint of the MF=10 grid at 150 MeV, an order of
magnitude in energy above the 14 MeV a D-T source produces and outside this
spectrum entirely. The converted `nuclide.arrow` files hash the same. Two
libraries carrying one evaluation give one answer, and anything that shows them
apart is reporting something other than the data.

### What else comes out of the same evaluations

`convert_to_arrow.py` writes four things, not one. The cross sections, a
`covariance.arrow` beside each of them, and two subsections under
`data/<source>/chain-<case>/`:

* `branching/`, which decides how much of an (n,2n) leaves the product in its
  metastable state rather than its ground state. That fraction is energy
  dependent, it lives in MF=8/9/10 of the same neutron evaluations, and for some
  foils it is the entire answer.
* `reactions/`, the network topology: which reaction on which parent gives which
  product, and the Q value of each. An evaluation states this outright, in the
  MT numbers it carries. TENDL-2025 gives Fe56 25 reaction channels where the
  endf-b8.1 chain gives 5, and W184 27 against 9.
* `covariance.arrow`, the MF=33 cross-section covariance, which is what puts the
  `+/- 6%` on the calculated decay heat. `--no-covariance` leaves it out; there
  is not much reason to, since it costs no measurable NJOY time and 4% more
  disk. An evaluation that carries no MF=33 writes none, and step 1 says which.

  It has to come from the ENDF, so a directory converted with `--no-covariance`
  cannot be given it without rerunning NJOY. Step 2 against one of those does
  not fail and does not report zero: it says no covariance was available and
  names the reconversion, and the report leaves both columns empty.

The first two matter for the same reason: every part of the answer that comes from
somewhere other than the library under test makes the C/E harder to read. Mean
deviation against the measurement on `2000exp_5min`, moving one subsection at a
time off endf-b8.1 and onto the library the cross sections came from:

| case | chain from endf-b8.1 | + branching | + topology | dominant product |
|---|---|---|---|---|
| Ti | 2.4% | 2.2% | **2.2%** | Sc50 |
| Fe | 5.3% | 5.7% | **6.0%** | Mn56 |
| Ag | 69.7% | 9.0% | **9.0%** | Ag108 |
| Nb | 81.4% | 10.8% | **10.8%** | **Nb94m**, 92% of the heat |
| W | 75.9% | 121.5% | **122.0%** | **W185m**, 97% of the heat |

Note which way some of these move. Fe gets worse at every step, and the last
column is the honest number for TENDL-2025 on iron precisely because it is the
one with the least endf-b8.1 left in it. A blend that happens to agree better
with the measurement is not evidence about either library.

The branching does the heavy lifting. Nb goes from 6x low to agreeing, Ag from
70% out to 9%, and the ground-state foils barely move, which is the tell: it
only matters where an isomer carries the heat. Read the first column as an
artefact and not as an endf-b8.1 result, because it pairs one library's
branching with another's cross sections.

The topology is worth much less on these five. The channels it adds are mostly
minor routes to products the foil already makes: iron's extra 20 channels move
the heat by 0.3%, all of it more Mn56 and Mn57 arriving by (n,d) and (n,np) on
Fe57 rather than only by (n,p) on Fe56, and tungsten's move it by up to 0.8%,
mostly more Ta185. Ti, Ag and Nb do not move at the printed precision. It is
worth doing anyway because it is one more thing the answer no longer borrows,
and on a foil whose heat runs through an unusual channel it would not be small.

Both subsections are built from the foil's own isotopes, which covers the
reactions that matter for a single 5 minute irradiation but is not the whole
chain. Both need decay data, branching to map a product's nuclear level to its
metastable state and the topology to know which nuclides exist, so the first run
fetches the ENDF/B-8.1 decay sublibrary (10.6 MB, 68 MB unpacked) into
`data/decay/`. Point `--decay-dir` at a copy you already have, or pass
`--no-branching` and `--no-reactions` to skip either and walk the columns back.

## The experiment

JAEA's Fusion Neutronics Source, as distributed through the IAEA CoNDERC decay
heat benchmark. A 1 g foil is held in a 14 MeV field of about 1.1e10 n/cm2/s,
then withdrawn, and its heat output is measured repeatedly as it cools.

### fns_data.json

The measured data in `fns_data.json` is reformatted from the benchmark archive
the IAEA publishes at [CoNDERC](https://nds.iaea.org/conderc/fusion/), which
`build_fns_data.py` reads directly:

```bash
python build_fns_data.py                    # fetches fns.zip, 14 MB, cached
python build_fns_data.py --fns-zip fns.zip  # a copy you already have
```

Nobody needs the archive to use this example. The script is kept because it is
the record of how the JSON was derived, and because re-deriving it is the only
way to check that record.

Three things go in, per experiment:

* the **spectrum**, 709 CCFE groups. The source files are written highest energy
  first, and are flipped once at build time to match the `CCFE-709` boundaries.
  There are only **four distinct spectra** across all 132 experiments, one per
  irradiation position in the FNS assembly, so the JSON stores them under
  `spectra_by_position` and each case names the position it sat at. Position 3
  covers all 73 foils of the 2000 campaign and 18 of the 1996 one.
* the **measurement**: time after shutdown, specific decay heat in microwatts
  per gram, and its uncertainty. Some source files are padded to a fixed length
  with all-zero rows, which are not measurements and are dropped. The times are
  cumulative from the end of the irradiation, in a unit the file never states,
  so it is inferred from the deck's cooling schedule and stored explicitly:
  minutes for the 5-minute experiments, days for the 7-hour one.
* the **FISPACT input deck**, for the five keywords that describe the physics:
  `DENSITY`, `MASS` and the element lines that follow it, `FLUX`, and `TIME`.
  The composition is the important one, and the only place it exists. The rest
  of a deck is solver settings and output requests, which mean nothing here, so
  no FISPACT input reader is needed.

The cooling steps yani takes are the intervals between measured points rather
than the deck's own steps, which is why the deck's cooling schedule is read at
build time but not carried into the JSON. The two are meant to agree and nearly
always do, but a few decks drift from their measurement (`Co/1996exp_5min` stops
at 54.7 minutes against a last measured point at 57.0). Stepping to the measured
times makes the comparison exact at every point for all 132 experiments.

## Which library supplies what

TENDL is a neutron-reaction library. Its sublibraries are neutron, proton,
deuteron, triton, He3, alpha, gamma, fission yields and thermal scattering.
There is no decay sublibrary, so it cannot say how long Mn56 lives or how much
energy each decay releases.

* From your library, converted here: the **reaction cross sections**, which set
  how fast each product is made, the **isomeric branching**, which sets which
  state it is made in, and the **reaction topology**, which says which reaction
  on which parent gives which product at all. All three come out of the same
  neutron evaluations, so all three move together when you change library.
* From endf-b8.1, downloaded by yani on first use: the **decay data**, meaning
  half-lives, decay modes and **decay energies**.

The split is not arbitrary, and it is drawn where the data runs out rather than
where it was convenient. The first group is everything a neutron evaluation
states and the second is what no neutron evaluation can, so this is the most of
a library swap that TENDL is able to express.

Fission yields are the one thing TENDL does publish that is not used here: none
of the 73 FNS foils is fissionable, so the subsection is never built.

### So what is a C/E here a statement about

Decay heat at time t is a sum over products of N(t) x lambda x Q. The library
under test sets N, how much of each product got made; endf-b8.1 sets lambda and
Q, how fast it decays and how much energy that releases. A C/E is therefore a
statement about production rates, measured against fixed decay data.

That is the useful decomposition for judging a neutron library, and it is worth
being explicit about because the measured quantity is decay heat, so the decay
energies scale the answer directly. Two consequences follow:

* A foil is only as good a test as its decay data is settled. Where the heat
  runs through one well known product, as iron's does through Mn56, the C/E is
  a clean read on one cross section. Where it does not, some of the deviation
  belongs to endf-b8.1.
* Changing library changes only N, so two libraries run through this example
  are compared like for like. W is the example worth looking at: TENDL-2025
  gives 122% and TENDL-2017 gives 84.5% on the same foil, against the same
  decay data, and 97% of that heat is W185m. That is a statement about how much
  W185m each library makes, and nothing else.
