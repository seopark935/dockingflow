"""Input parsing and validation for dockingflow.

This module has no side effects beyond reading files: it turns
`setup.txt`, `tranches.txt`, `recList.txt`, and `geoList.txt` into typed,
validated Python objects that `pipeline.py` (and `gui.py`) consume. Keeping
all format-specific parsing and validation logic here means both the CLI
and the GUI fail on the same bad input in the same way, with the same
error messages.

File formats
------------
- setup.txt: one `KEY=VALUE` pair per line; `#`-prefixed lines and blank
  lines are ignored. See `load_setup`.
- tranches.txt (the "map" file): a TSV with header
  `curl_script  log_p  molecular_weight`, one row per ZINC tranche. See
  `load_tranches_tsv`.
- recList.txt: one receptor `.pdbqt` path per line, resolved relative to
  the file's own directory.
- geoList.txt: one grid box per line, `center_x center_y center_z size_x
  size_y size_z`, paired line-for-line with recList.txt. See
  `load_docking_targets`.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

@dataclass(frozen=True)
class Tranche:
    """One row of tranches.txt: a ZINC tranche identified by (log_p, molecular_weight)."""

    curl_script: Path
    log_p: float
    molecular_weight: int

@dataclass(frozen=True)
class GridBox:
    """An AutoDock Vina-style search box: a center point and box dimensions, in angstroms."""

    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float

@dataclass(frozen=True)
class DockingTarget:
    """One receptor paired with the grid box docking should search within it."""

    receptor: Path
    box: GridBox

def load_setup(setup_path: Path) -> Dict[str, str]:
    """Parse a `KEY=VALUE`-per-line setup file into a plain string dict.

    Blank lines and lines starting with `#` are skipped. Every other line
    must contain exactly one `=` with non-empty key and value (values may
    themselves contain `=`, since only the first `=` is used as the split
    point).

    Args:
        setup_path: Path to the setup file (typically `setup.txt`).

    Returns:
        A dict of raw string values — callers are responsible for parsing
        numeric fields themselves (see `validate_inputs`).

    Raises:
        FileNotFoundError: If `setup_path` doesn't exist.
        ValueError: If any non-blank, non-comment line isn't valid
            `KEY=VALUE`.
    """
    setup_path = setup_path.expanduser().resolve()
    if not setup_path.exists():
        raise FileNotFoundError(f"setup file not found: {setup_path}")

    setup: Dict[str, str] = {}
    for raw in setup_path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            raise ValueError(f"Bad setup line (expected KEY=VALUE): {raw}")
        key, val = line.split("=", 1)
        key = key.strip()
        val = val.strip()
        if not key or not val:
            raise ValueError(f"Bad setup line (empty key/value): {raw}")
        setup[key] = val

    return setup

def load_tranches_tsv(map_path: Path) -> List[Tranche]:
    """Parse the tranche mapping file into a list of `Tranche` objects.

    Expects a header row of exactly `curl_script log_p molecular_weight`
    (whitespace-separated) followed by one data row per tranche. Each
    `curl_script` entry is resolved relative to the mapping file's own
    directory, not the current working directory, so the mapping file can
    be run from anywhere.

    Args:
        map_path: Path to the mapping file (typically `tranches.txt`).

    Returns:
        One `Tranche` per data row, in file order.

    Raises:
        FileNotFoundError: If `map_path` doesn't exist.
        ValueError: If the file is empty, the header doesn't match exactly,
            a data row doesn't have exactly 3 fields, or `log_p`/
            `molecular_weight` aren't parseable as float/int respectively.
    """
    map_path = map_path.expanduser().resolve()
    if not map_path.exists():
        raise FileNotFoundError(f"mapping file not found: {map_path}")

    lines = [ln.rstrip("\n") for ln in map_path.read_text().splitlines() if ln.strip()]
    if not lines:
        raise ValueError("mapping file is empty")

    header = lines[0].split()
    expected = ["curl_script", "log_p", "molecular_weight"]
    if header != expected:
        raise ValueError(f"Bad header. Expected {expected} but got {header}")

    tranches: List[Tranche] = []
    for i, ln in enumerate(lines[1:], start=2):
        parts = ln.split()
        if len(parts) != 3:
            raise ValueError(f"Line {i}: expected 3 tab-separated fields, got {len(parts)}: {ln}")
        curl_script_s, logp_s, mw_s = parts
        try:
            logp = float(logp_s)
        except ValueError:
            raise ValueError(f"Line {i}: log_p is not a float: {logp_s}")
        try:
            mw = int(mw_s)
        except ValueError:
            raise ValueError(f"Line {i}: molecular_weight is not an int: {mw_s}")

        curl_script = (map_path.parent / curl_script_s).expanduser().resolve()
        tranches.append(Tranche(curl_script=curl_script, log_p=logp, molecular_weight=mw))

    return tranches

def load_docking_targets(setup: Dict[str, str]) -> List[DockingTarget]:
    """Pair up recList.txt lines (receptor paths) with geoList.txt lines (grid boxes).

    Both files are line-oriented lists so that multiple receptors/binding
    sites can be docked against in a single run; line i of recList is
    paired with line i of geoList. Receptor paths are resolved relative to
    recList.txt's own directory.

    Args:
        setup: The parsed setup dict, as returned by `load_setup`. Must
            contain `recList` and `geoList` keys pointing at existing files.

    Returns:
        One `DockingTarget` per line pair, in file order.

    Raises:
        ValueError: If recList is empty, if recList and geoList have
            different numbers of non-blank lines, or if a geoList line
            doesn't parse as exactly 6 floats.
        FileNotFoundError: If a receptor path listed in recList doesn't
            exist on disk.
    """
    rec_path = Path(setup["recList"]).expanduser().resolve()
    geo_path = Path(setup["geoList"]).expanduser().resolve()

    rec_lines = [ln.strip() for ln in rec_path.read_text().splitlines() if ln.strip()]
    geo_lines = [ln.strip() for ln in geo_path.read_text().splitlines() if ln.strip()]

    if not rec_lines:
        raise ValueError(f"recList is empty: {rec_path}")
    if len(rec_lines) != len(geo_lines):
        raise ValueError(
            f"recList ({len(rec_lines)} entries) and geoList ({len(geo_lines)} entries) "
            "must have the same number of lines"
        )

    targets: List[DockingTarget] = []
    for i, (rec_line, geo_line) in enumerate(zip(rec_lines, geo_lines), start=1):
        receptor = (rec_path.parent / rec_line).expanduser().resolve()
        if not receptor.exists():
            raise FileNotFoundError(f"recList line {i}: receptor not found: {receptor}")

        parts = geo_line.split()
        if len(parts) != 6:
            raise ValueError(
                f"geoList line {i}: expected 6 numbers (cx cy cz sx sy sz), got {len(parts)}: {geo_line}"
            )
        try:
            cx, cy, cz, sx, sy, sz = (float(p) for p in parts)
        except ValueError:
            raise ValueError(f"geoList line {i}: could not parse floats: {geo_line}")

        targets.append(DockingTarget(receptor=receptor, box=GridBox(cx, cy, cz, sx, sy, sz)))

    return targets

def validate_inputs(setup: Dict[str, str], tranches: List[Tranche], workdir: Path) -> None:
    """Fail fast on any malformed input before the pipeline does real work.

    Checks, in order: required setup.txt keys are present; recList/geoList/
    ligList point at existing files; filter_percent is a number in (0, 100];
    cores is an integer >= 1; energy_range (if present) is a positive
    number; recList/geoList pair up correctly (via `load_docking_targets`);
    at least one tranche was loaded; every tranche's curl_script exists and
    is unique; no two tranches share a (log_p, molecular_weight) identity;
    and workdir exists.

    Note that `ligList` is only checked for existence here — its contents
    aren't consumed directly by the pipeline. Per-tranche ligand lists are
    generated automatically during the unpack stage
    (`pipeline.unpack_tranche`) from each tranche's downloaded ligands.

    Args:
        setup: Parsed setup dict, as returned by `load_setup`.
        tranches: Parsed tranche list, as returned by `load_tranches_tsv`.
        workdir: The run's output directory; must already exist.

    Raises:
        ValueError: For any structurally invalid input (see checks above).
        FileNotFoundError: For any referenced path that doesn't exist.
    """
    required = ["recList", "geoList", "ligList", "filter_percent", "cores"]
    missing = [k for k in required if k not in setup]
    if missing:
        raise ValueError(f"setup missing required keys: {missing}")

    # Validate paths
    rec = Path(setup["recList"]).expanduser().resolve()
    geo = Path(setup["geoList"]).expanduser().resolve()
    lig = Path(setup["ligList"]).expanduser().resolve()

    for p, name in [(rec, "recList"), (geo, "geoList"), (lig, "ligList")]:
        if not p.exists():
            raise FileNotFoundError(f"{name} does not exist: {p}")
        if not p.is_file():
            raise ValueError(f"{name} is not a file: {p}")

    # Validate numeric parameters
    try:
        fp = float(setup["filter_percent"])
    except ValueError:
        raise ValueError("filter_percent must be a number")
    if not (0.0 < fp <= 100.0):
        raise ValueError("filter_percent must be in (0, 100]")

    try:
        cores = int(setup["cores"])
    except ValueError:
        raise ValueError("cores must be an integer")
    if cores < 1:
        raise ValueError("cores must be >= 1")

    if "energy_range" in setup:
        try:
            er = float(setup["energy_range"])
        except ValueError:
            raise ValueError("energy_range must be a number")
        if er <= 0:
            raise ValueError("energy_range must be > 0")

    # Fail fast on malformed recList/geoList pairing rather than during docking.
    load_docking_targets(setup)

    if not tranches:
        raise ValueError("No tranches loaded from mapping file")

    # Validate tranche scripts + uniqueness
    seen_pairs: set[Tuple[float, int]] = set()
    seen_scripts: set[Path] = set()
    for t in tranches:
        if t.curl_script in seen_scripts:
            raise ValueError(f"Duplicate curl_script in mapping: {t.curl_script}")
        seen_scripts.add(t.curl_script)

        if not t.curl_script.exists():
            raise FileNotFoundError(f"curl script not found: {t.curl_script}")

        key = (t.log_p, t.molecular_weight)
        if key in seen_pairs:
            raise ValueError(f"Duplicate tranche (log_p, molecular_weight) in mapping: {key}")
        seen_pairs.add(key)

    # Workdir sanity
    if not workdir.exists():
        raise FileNotFoundError(f"workdir does not exist: {workdir}")

def format_summary(setup: Dict[str, str], tranches: List[Tranche], workdir: Path) -> str:
    """Render a short human-readable summary of a validated run configuration.

    Shown to the user right before Stage 0 starts (CLI) or as the result of
    clicking "Validate" (GUI), so they can sanity-check the run before any
    files are downloaded.
    """
    logps = sorted({t.log_p for t in tranches})
    mws = sorted({t.molecular_weight for t in tranches})
    return (
        "Inputs look valid.\n"
        f"Workdir: {workdir}\n"
        f"Tranches: {len(tranches)}\n"
        f"log_p bins ({len(logps)}): {logps}\n"
        f"MW bins ({len(mws)}): {mws}\n"
        f"filter_percent: {setup['filter_percent']}\n"
        f"cores: {setup['cores']}\n"
        f"recList: {setup['recList']}\n"
        f"geoList: {setup['geoList']}\n"
        f"ligList: {setup['ligList']}\n"
    )
