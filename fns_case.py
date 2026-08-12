"""Reading one FNS decay-heat case out of `fns_data.json`.

That file holds the whole benchmark: 73 foils, 132 experiments, and the four
distinct neutron spectra they were measured in. `build_fns_data.py` packs it
from the benchmark's own files and says how.

A case is named for its foil, which is an element symbol for most of them and
an alloy like SS316 for a few, plus which experiment to take.
"""

import functools
import json
import pathlib

import numpy as np

GROUPS = 709

# Which experiment to run when none is asked for. The 2000 campaign measured
# all 73 foils; the two 1996 ones cover about 30 each, so this is the only
# choice that works for every case. The alternatives are worth knowing about:
# 1996exp_5min has slightly tighter uncertainties (5.0% against 6.8% median),
# and 1996exp_7hour irradiates for 7.6 hours and follows the foil out to 400
# days, which tests long-lived products the 5-minute cases never reach.
PREFERRED = ["2000exp_5min", "1996exp_5min", "1996exp_7hour"]


def default_path():
    return pathlib.Path(__file__).resolve().parent / "fns_data.json"


@functools.lru_cache(maxsize=4)
def _document(path):
    if not path.is_file():
        raise SystemExit(f"{path} is missing. Rebuild it with build_fns_data.py.")
    return json.loads(path.read_text())


class Case:
    """One foil, one experiment: what to irradiate, how, and what was measured."""

    def __init__(self, name, experiment, entry, spectrum):
        self.name = name
        self.experiment = experiment
        self.density = entry["density_g_per_cm3"]
        self.mass_g = entry["mass_g"]
        self.composition = entry["composition_mass_fraction"]
        self.irradiation = [(seconds, flux) for seconds, flux in entry["irradiation_s_and_flux"]]
        self.position = entry["position"]
        self.time_unit = entry["time_unit"]
        self.spectrum = np.array(spectrum)
        self.times = np.array(entry["times"])
        self.measured = np.array(entry["measured_uW_per_g"])
        self.uncertainty = np.array(entry["uncertainty_uW_per_g"])

        # The deck's own cooling steps are not carried: yani cools in the
        # intervals between measured points instead, so the comparison is exact
        # at every point. A few decks drift from their own measurement
        # (Co/1996exp_5min stops at 54.7 minutes against a last point at 57.0).
        seconds = {"minutes": 60.0, "hours": 3600.0, "days": 86400.0}[self.time_unit]
        self.cooling = np.diff(np.concatenate(([0.0], self.times * seconds))).tolist()
        if min(self.cooling) <= 0.0:
            raise SystemExit(f"{name}/{experiment}: measured times are not increasing")

    @property
    def total_flux(self):
        """Total flux [n/cm2/s], which the irradiation pulse carries."""
        return self.irradiation[0][1]

    @property
    def elements(self):
        return sorted(self.composition)

    def describe(self):
        parts = ", ".join(
            f"{el} {frac * 100:.4g}%" for el, frac in sorted(self.composition.items())
        )
        return (f"{self.name} / {self.experiment}: {self.mass_g:g} g at "
                f"{self.density:g} g/cm3 ({parts}), "
                f"{sum(t for t, _ in self.irradiation):g} s at "
                f"{self.total_flux:.4g} n/cm2/s at position {self.position}, "
                f"{len(self.cooling)} cooling steps")


def cases(path=None):
    """Every foil in the benchmark."""
    return sorted(_document(path or default_path())["cases"])


def experiments(case, path=None):
    """Experiment names available for one foil."""
    return sorted(_document(path or default_path())["cases"].get(case, {}))


def load(case, experiment=None, path=None):
    """Assemble one Case."""
    document = _document(path or default_path())
    available = document["cases"].get(case)
    if not available:
        raise SystemExit(
            f"no foil {case!r} in the benchmark. Available: {', '.join(cases(path))}"
        )
    if experiment is None:
        experiment = next((e for e in PREFERRED if e in available), sorted(available)[0])
    elif experiment not in available:
        raise SystemExit(
            f"{case} has no experiment {experiment!r}. "
            f"Available: {', '.join(sorted(available))}"
        )
    entry = available[experiment]
    return Case(case, experiment, entry, document["spectra_by_position"][entry["position"]])
