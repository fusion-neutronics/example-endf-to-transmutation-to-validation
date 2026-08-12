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
   plotted against the measurement.

It defaults to iron, and `--case` takes any of the 73 FNS foils.

There is no neutron transport anywhere in this example. The FNS experiment
published the spectrum the foil actually saw, so that spectrum is the input and
yani collapses it against the cross sections to get one-group reaction rates.
That is what yani is for: transmutation and activation given a spectrum.

## Install

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

That is the whole install, and it no longer needs a checkout of anything.
`yani` both reads the nuclear data and makes it: the ENDF conversion, the
reaction topology and the isomeric branching are all functions on the wheel.

It installs from git rather than PyPI because no wheel is published under this
name yet. The package was called `yats` until recently, and that name is
unusable on PyPI: it belongs to an unrelated Twitter scraper, so
`pip install yats` succeeds and gives you a package with none of these
functions. `yani` is free, so this becomes a plain `pip install` once the first
wheel is up.

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
are written to disk. Conversion runs about twenty seconds per isotope and the
results are cached, so the second run of the script does nothing.

`run_transmutation.py` writes `results/<source>/fns_<case>_<experiment>.png` and
a JSON of the same numbers, including the per-nuclide breakdown of the heat.

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

Swapping the cross sections for ENDF/B-8.1, with everything else held fixed,
turns this into a library comparison:

```bash
python convert_to_arrow.py --endf-dir /path/to/endfb-8.1 --library endf-b8.1 \
    --output data/b81
python run_transmutation.py --cross-sections data/b81 --output results/b81
```

### What else comes out of the same evaluations

`convert_to_arrow.py` writes three things, not one. The cross sections, and two
subsections under `data/chain/`:

* `branching/`, which decides how much of an (n,2n) leaves the product in its
  metastable state rather than its ground state. That fraction is energy
  dependent, it lives in MF=8/9/10 of the same neutron evaluations, and for some
  foils it is the entire answer.
* `reactions/`, the network topology: which reaction on which parent gives which
  product, and the Q value of each. An evaluation states this outright, in the
  MT numbers it carries. TENDL-2025 gives Fe56 25 reaction channels where the
  endf-b8.1 chain gives 5, and W184 27 against 9.

Both matter for the same reason: every part of the answer that comes from
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
