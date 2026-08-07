# ENDF to transmutation to validation

An end-to-end example: raw evaluated nuclear data goes in one end, and a
comparison against a measured experiment comes out the other.

The whole chain is three steps.

1. **ENDF in.** Evaluations for the natural isotopes of whatever the foil is
   made of. TENDL-2025 by default, any library you point it at otherwise.
2. **Arrow out.** NJOY reconstructs the resonances and Doppler broadens to
   294 K, and `nuclear_data_to_arrow` writes the result as Arrow, which is the
   form yats reads. This happens on your machine; nothing pre-processed is
   downloaded.
3. **Validation.** yats irradiates a 1 g foil with the measured neutron
   spectrum from the JAEA FNS decay-heat experiment, follows the activation
   products through the cooling schedule, and the specific decay heat is
   plotted against the measurement.

It defaults to iron, and `--case` takes any of the 73 FNS foils.

There is no neutron transport anywhere in this example. The FNS experiment
published the spectrum the foil actually saw, so that spectrum is the input and
yats collapses it against the cross sections to get one-group reaction rates.
That is what yats is for: transmutation and activation given a spectrum, with
transport left to yamc or to an experiment.

## Install

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```


* `yats` and `nuclear_data_to_arrow` come from the yamc checkout, expected at
  `../yamc-org/yamc`. The relative paths resolve against the directory pip runs
  in, so install from this one. yats compiles a Rust extension, which is the
  slow part and needs a Rust toolchain (https://rustup.rs); everything else is
  a wheel.
* the `njoy2016` wheel is on `https://shimwell.github.io/wheels`, not PyPI. It
  puts an `njoy` binary in the venv, which is what the conversion shells out to.
* the `endf` reader is the `local-develop` branch of
  `github.com/shimwell/endf-python`.

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

With no `--endf-dir` it fetches TENDL-2025 for you:

```bash
python convert_to_arrow.py --tarball /path/to/TENDL-n.tgz
python convert_to_arrow.py           # streams TENDL-n.tgz from tendl.imperial.ac.uk
```

The last form transfers 3.5 GB whatever you asked for, because TENDL publishes
its neutron sublibrary as a single archive; only the wanted members are written
to disk. Conversion runs about twenty seconds per isotope and the results are
cached, so the second run of the script does nothing.

`run_transmutation.py` writes `results/fns_<case>_<experiment>.png` and a JSON
of the same numbers, including the per-nuclide breakdown of the heat.

## What comes out

For iron on the default experiment, yats tracks the measurement to 5.3% mean
deviation, starting on top of it and drifting to 8% high by the end of the
hour. Run the same iron against `1996exp_5min` and it comes out 7.8% low
instead, C/E 0.886 rising to 0.956. The two campaigns bracket the calculation,
which is a fair reminder that a single C/E number is a statement about one
measurement and not only about the data.

The breakdown says why the shape is what it is. Mn56 from Fe56(n,p), half-life
2.58 hours, is 81% of the heat at the first point and 99.9% at the last. The
rest of the early heat is Mn57 and Fe53, both of which are gone within minutes,
so the first few points test three cross sections and the tail tests one.

Swapping the cross sections for ENDF/B-8.1, with everything else held fixed,
turns this into a library comparison:

```bash
python convert_to_arrow.py --endf-dir /path/to/endfb-8.1 --library endf-b8.1 \
    --output data/b81
python run_transmutation.py --cross-sections data/b81 --output results/b81
```

### Foils whose heat comes from an isomer

Accuracy is not uniform across the 73 foils, and the split is not random. On
`2000exp_5min`, against the same benchmark run through the full TENDL-2025
chain:

| case | this example | full TENDL-2025 chain | dominant product |
|---|---|---|---|
| Ti | 2.4% | 2.2% | Sc50 |
| Fe | 5.3% | 6.0% | Mn56 |
| Cu | 5.6% | 6.3% | Cu62 |
| Al | 7.1% | 7.1% | Al28 |
| Ag | 69.7% | 9.0% | Ag108 |
| Nb | 81.4% | 10.8% | **Nb94m**, 92% of the heat |
| W | 75.9% | 122.0% | **W185m**, 97% of the heat |

Every case that agrees is one whose decay heat comes from a ground-state
product. Every case that does not is one where an isomer carries almost all of
it, and Nb comes out 6x low.

That is this example's one real approximation showing itself. It takes the
**cross sections** from your library but holds the **transmutation network**,
including the isomeric branching ratios, on endf-b8.1. How much of an (n,2n)
lands in the metastable state rather than the ground state is energy dependent,
a scalar branching table cannot express it, and for these foils it decides the
answer. The full chain carries a `branching/` subsection extracted from the
same TENDL evaluations as the cross sections, which is what closes the gap.

Building that subsection here would mean pulling in the ENDF/B-8.1 decay
sublibrary as well, so it is left out. Read the isomer-dominated foils as a
demonstration of what the branching data is for, not as a yats result.

## The experiment

JAEA's Fusion Neutronics Source, as distributed through the IAEA CoNDERC decay
heat benchmark. A 1 g foil is held in a 14 MeV field of about 1.1e10 n/cm2/s,
then withdrawn, and its heat output is measured repeatedly as it cools.

### fns_data.json

The measured data is redistributed here, not generated. It arrives as 396 files
across 73 foil directories, three per experiment, via
[jbae11/openmc_activator](https://github.com/jbae11/openmc_activator), whose own
comparison notebook is worth reading alongside this one. `build_fns_data.py`
packs that tree into one 175 KB `fns_data.json`, which is what gets committed
and what `fns_case.py` reads:

```bash
python build_fns_data.py --fns-dir ~/openmc_activator/fns
```

Nobody needs the original tree to use this example. The script is kept because
it is the record of how the JSON was derived, and because re-deriving it is the
only way to check that record.

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

The cooling steps yats takes are the intervals between measured points rather
than the deck's own steps, which is why the deck's cooling schedule is read at
build time but not carried into the JSON. The two are meant to agree and nearly
always do, but a few decks drift from their measurement (`Co/1996exp_5min` stops
at 54.7 minutes against a last measured point at 57.0). Stepping to the measured
times makes the comparison exact at every point for all 132 experiments.

## Which library supplies what

TENDL is a neutron-reaction library. It has no decay sublibrary, so it cannot
say how long Mn56 lives or how much energy each decay releases, and it cannot
on its own define the transmutation network.

* From TENDL-2025, converted here: the **reaction cross sections**, which set
  how fast each product is made.
* From endf-b8.1, downloaded by yats on first use: the **decay data**, the
  **decay energies**, and the **reaction network** that says which reaction on
  which parent gives which product.

Swapping `--cross-sections` for a directory built from another library, with
everything else held fixed, turns this example into a library comparison.
