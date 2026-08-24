"""Which nuclear data a run used, as one folder name.

Every artefact a run produces is kept under that name: the unpacked ENDF, the
converted Arrow, and the results. Two libraries therefore sit side by side
rather than overwriting each other, which is what makes comparing them possible
at all.

The collision this avoids is not hypothetical. TENDL names its evaluations the
same way in every release, so ``n-Fe056.tendl`` from TENDL-2017 and from
TENDL-2025 have the same basename. Unpacked into one directory the second
library silently reuses the first one's files, and a "TENDL-2017" sweep quietly
reports TENDL-2025 numbers.
"""

from __future__ import annotations

import pathlib
import re

# The neutron sublibrary tarballs we know how to fetch. TENDL publishes one
# archive per release at a stable path, so the release year is the only thing
# that changes. Both are around 3 GB.
#
# The two archives are laid out differently: TENDL-2025 keeps its evaluations
# flat at the archive root (``n-Fe056.tendl``), TENDL-2017 nests them
# (``neutron_file/Fe/Fe056/lib/endf/n-Fe056.tendl``). The extractor matches on
# the basename, so both work without special-casing.
TENDL_TARBALLS = {
    "tendl-2017": "https://tendl.imperial.ac.uk/tendl_2017/tar_files/TENDL-n.tgz",
    "tendl-2025": "https://tendl.imperial.ac.uk/tendl_2025/tar_files/TENDL-n.tgz",
}

DEFAULT_LIBRARY = "tendl-2025"


def _sanitise(name: str) -> str:
    """A string safe to use as one path component."""
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "-", name).strip("-.")
    return cleaned or "endf"


def slug(library: str, endf_dir=None, override=None) -> str:
    """Folder name identifying the nuclear data a run used.

    A downloaded library is named for itself (``tendl-2017``). Evaluations of
    your own are named for the directory they came from, so two directories of
    your own stay separate too, while ``--library`` still reaches the Arrow
    provenance as the label the data is stamped with.

    `override` (``--source``) wins over both, for when the directory name is not
    what you want the results filed under: a directory called
    ``tendl-2025-endf`` holds exactly the library ``tendl-2025``, and there is
    no reason for the results to say otherwise.
    """
    if override:
        return _sanitise(override)
    if endf_dir is not None:
        return _sanitise(pathlib.Path(endf_dir).expanduser().resolve().name)
    return _sanitise(library)


def tarball_url(library: str) -> str:
    """Where to fetch `library` from, or a message naming the alternatives."""
    try:
        return TENDL_TARBALLS[library]
    except KeyError:
        known = ", ".join(sorted(TENDL_TARBALLS))
        raise SystemExit(
            f"no download is known for library {library!r}. Either pass one of "
            f"{known}, or point --endf-dir at evaluations you already have and "
            f"keep --library {library!r} as the label they are stamped with."
        ) from None
