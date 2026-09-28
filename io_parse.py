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
  lines are ignored. Relative paths are resolved against setup.txt's own
  directory. See `load_setup`.
- tranches.txt (the "map" file): a whitespace-separated table, one row per
  ZINC tranche, with either header `curl_script  tranche` (ZINC22 tranche
  codes like `H04M000`, as written by `zinc_split.py`) or the older
  `curl_script  log_p  molecular_weight` (ZINC20-style). See
  `load_tranches_tsv`.
- recList.txt: one receptor `.pdbqt` path per line, resolved relative to
  the file's own directory.
- geoList.txt: one grid box per line, `center_x center_y center_z size_x
  size_y size_z`, paired line-for-line with recList.txt. See
  `load_docking_targets`.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ZINC22 3D tranche code: H<heavy atom count><M|P><logP code>, e.g. H04M000,
# H17P050. The logP code is logP * 100 (M = minus, P = plus).
ZINC22_TRANCHE_RE = re.compile(r"H(\d{2})([MP])(\d{3})")

# setup.txt keys holding file paths, resolved relative to setup.txt itself.
SETUP_PATH_KEYS = ("recList", "geoList", "ligList")


@dataclass(frozen=True)
class Tranche:
    """One row of tranches.txt.

    A tranche is identified either by a ZINC22 tranche `code` (with
    `heavy_atoms`/`log_p` decoded from it) or, for the older map format, by
    its `(log_p, molecular_weight)` pair.
    """

    curl_script: Path
    log_p: float
    molecular_weight: Optional[int] = None
    heavy_atoms: Optional[int] = None
    code: Optional[str] = None


def parse_zinc22_code(code: str) -> Tuple[int, float]:
    """Decode a ZINC22 tranche code like `H04M000` into `(heavy_atoms, log_p)`.

    Raises:
        ValueError: If `code` isn't a ZINC22 tranche code.
    """
    m = ZINC22_TRANCHE_RE.fullmatch(code)
    if not m:
        raise ValueError(f"Not a ZINC22 tranche code (expected e.g. H04M000): {code}")
    hac, sign, logp = m.groups()
    return int(hac), (-1 if sign == "M" else 1) * int(logp) / 100

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

    Path-valued keys (`recList`, `geoList`, `ligList`) are resolved relative
    to setup.txt's own directory and returned as absolute paths, so a run
    behaves the same no matter which directory it's launched from.

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

    for key in SETUP_PATH_KEYS:
        if key in setup:
            setup[key] = str((setup_path.parent / Path(setup[key]).expanduser()).resolve())

    return setup

def load_tranches_tsv(map_path: Path) -> List[Tranche]:
    """Parse the tranche mapping file into a list of `Tranche` objects.

    Expects a header row (whitespace-separated) of exactly either:

      - `curl_script tranche`: ZINC22 tranche codes such as `H04M000`, as
        written by `zinc_split.py`; or
      - `curl_script log_p molecular_weight`: the older ZINC20-style format.

    followed by one data row per tranche. Each `curl_script` entry is
    resolved relative to the mapping file's own directory, not the current
    working directory, so the mapping file can be run from anywhere.

    Args:
        map_path: Path to the mapping file (typically `tranches.txt`).

    Returns:
        One `Tranche` per data row, in file order.

    Raises:
        FileNotFoundError: If `map_path` doesn't exist.
        ValueError: If the file is empty, the header matches neither format,
            a data row has the wrong number of fields, or a field doesn't
            parse (bad tranche code, or `log_p`/`molecular_weight` not
            float/int).
    """
    map_path = map_path.expanduser().resolve()
    if not map_path.exists():
        raise FileNotFoundError(f"mapping file not found: {map_path}")

    lines = [ln.rstrip("\n") for ln in map_path.read_text().splitlines() if ln.strip()]
    if not lines:
        raise ValueError("mapping file is empty")

    header = lines[0].split()
    zinc22_header = ["curl_script", "tranche"]
    legacy_header = ["curl_script", "log_p", "molecular_weight"]
    if header not in (zinc22_header, legacy_header):
        raise ValueError(f"Bad header. Expected {zinc22_header} or {legacy_header} but got {header}")

    tranches: List[Tranche] = []
    if header == zinc22_header:
        for i, ln in enumerate(lines[1:], start=2):
            parts = ln.split()
            if len(parts) != 2:
                raise ValueError(f"Line {i}: expected 2 fields (curl_script tranche), got {len(parts)}: {ln}")
            curl_script_s, code = parts
            try:
                hac, logp = parse_zinc22_code(code)
            except ValueError as exc:
                raise ValueError(f"Line {i}: {exc}")
            curl_script = (map_path.parent / curl_script_s).expanduser().resolve()
            tranches.append(Tranche(curl_script=curl_script, log_p=logp, heavy_atoms=hac, code=code))
        return tranches

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

# VinaLC's defaults, which setup.txt can override.
DEFAULT_EXHAUSTIVENESS = 8
DEFAULT_GRANULARITY = 0.375

# Memory model for one VinaLC worker rank. Each docking job builds a Vina
# grid map (double precision) for every atom type in the ligand, spanning
# the whole search box at `granularity` spacing; drug-like ligands use up to
# ~8 of Vina's XS atom types. On top of that is a fixed per-process overhead
# (receptor, precalculated tables, MPI buffers). The master rank only hands
# out work and is covered by the same overhead figure.
GRID_ATOM_TYPES = 8
RANK_OVERHEAD_GB = 0.2


def estimate_rank_memory_gb(targets: List["DockingTarget"], granularity: float = DEFAULT_GRANULARITY) -> float:
    """Estimated peak memory of one VinaLC worker, in GB, for the largest grid box in `targets`."""
    points = max(
        math.prod(math.ceil(size / granularity) + 1 for size in (t.box.size_x, t.box.size_y, t.box.size_z))
        for t in targets
    )
    return points * 8 * GRID_ATOM_TYPES / 1e9 + RANK_OVERHEAD_GB


def plan_mpi_ranks(cores: int, exhaustiveness: int, memory_gb: Optional[float], rank_memory_gb: float) -> int:
    """How many MPI ranks (1 master + workers) fit in `cores` and `memory_gb`.

    Each VinaLC worker runs ~`exhaustiveness` search threads, so a core
    budget supports `cores // exhaustiveness` workers; a memory budget
    supports `(memory_gb - master) // rank_memory_gb`. Always at least one
    worker, since VinaLC refuses to run with fewer than 2 ranks.
    """
    workers = cores // exhaustiveness
    if memory_gb is not None:
        workers = min(workers, int((memory_gb - RANK_OVERHEAD_GB) // rank_memory_gb))
    return max(1, workers) + 1


@dataclass(frozen=True)
class VinaLCOptions:
    """Docking settings from setup.txt, as passed to the VinaLC command line.

    `mpi_ranks` is the `-np` given to mpirun. VinaLC runs one master rank
    (which only hands out work) plus workers, and each worker runs roughly
    `exhaustiveness` search threads of its own — so launching one rank per
    core oversubscribes the machine ~`exhaustiveness`-fold. Unless setup.txt
    sets `mpi_ranks` explicitly, it's derived from the `cores` and optional
    `memory_gb` budgets by `plan_mpi_ranks`.
    """

    mpi_ranks: int
    energy_range: str
    exhaustiveness: int
    num_modes: int
    granularity: Optional[str] = None
    seed: Optional[str] = None

    def cli_args(self) -> List[str]:
        """VinaLC flags for these options (the list-file flags are added separately)."""
        args = [
            "--exhaustiveness", str(self.exhaustiveness),
            "--num_modes", str(self.num_modes),
            "--energy_range", self.energy_range,
        ]
        if self.granularity is not None:
            args += ["--granularity", self.granularity]
        if self.seed is not None:
            args += ["--seed", self.seed]
        return args


def _positive_number(setup: Dict[str, str], key: str) -> float:
    try:
        val = float(setup[key])
    except ValueError:
        raise ValueError(f"{key} must be a number")
    if val <= 0:
        raise ValueError(f"{key} must be > 0")
    return val


def _positive_int(setup: Dict[str, str], key: str, minimum: int = 1) -> int:
    try:
        val = int(setup[key])
    except ValueError:
        raise ValueError(f"{key} must be an integer")
    if val < minimum:
        raise ValueError(f"{key} must be >= {minimum}")
    return val


def load_vinalc_options(setup: Dict[str, str]) -> VinaLCOptions:
    """Parse and validate setup.txt's docking keys into a `VinaLCOptions`.

    Keys: `cores` (required), and optionally `memory_gb`, `mpi_ranks`,
    `exhaustiveness` (default 8), `num_modes` (default 9), `energy_range`
    (default 3), `granularity`, `seed`. When `memory_gb` is set, the grid
    boxes from recList/geoList are read to estimate per-rank memory.

    Raises:
        ValueError: If any of those keys is present but malformed.
    """
    cores = _positive_int(setup, "cores")
    exhaustiveness = _positive_int(setup, "exhaustiveness") if "exhaustiveness" in setup else DEFAULT_EXHAUSTIVENESS
    num_modes = _positive_int(setup, "num_modes") if "num_modes" in setup else 9

    energy_range = setup.get("energy_range", "3")
    _positive_number({"energy_range": energy_range}, "energy_range")

    granularity = setup.get("granularity")
    if granularity is not None:
        _positive_number(setup, "granularity")

    memory_gb = _positive_number(setup, "memory_gb") if "memory_gb" in setup else None

    if "mpi_ranks" in setup:
        # VinaLC refuses to run with fewer than 2 processes (master + 1 worker).
        mpi_ranks = _positive_int(setup, "mpi_ranks", minimum=2)
    else:
        rank_gb = 0.0
        if memory_gb is not None:
            rank_gb = estimate_rank_memory_gb(
                load_docking_targets(setup), float(granularity or DEFAULT_GRANULARITY)
            )
        mpi_ranks = plan_mpi_ranks(cores, exhaustiveness, memory_gb, rank_gb)

    seed = setup.get("seed")
    if seed is not None:
        try:
            int(seed)
        except ValueError:
            raise ValueError("seed must be an integer")

    return VinaLCOptions(
        mpi_ranks=mpi_ranks,
        energy_range=energy_range,
        exhaustiveness=exhaustiveness,
        num_modes=num_modes,
        granularity=granularity,
        seed=seed,
    )


def tranche_label(t: Tranche) -> str:
    """Filesystem-safe, deterministic tranche label (the ZINC22 code when there is one)."""
    if t.code:
        return t.code
    logp = f"{t.log_p:.2f}".rstrip("0").rstrip(".")
    return f"LP{logp}_MW{t.molecular_weight}"


def tranche_size(t: Tranche) -> str:
    """Human-readable size bin: heavy atom count (ZINC22) or molecular weight (legacy)."""
    return f"HAC {t.heavy_atoms}" if t.heavy_atoms is not None else f"MW {t.molecular_weight}"


def validate_inputs(setup: Dict[str, str], tranches: List[Tranche], workdir: Path) -> None:
    """Fail fast on any malformed input before the pipeline does real work.

    Checks, in order: required setup.txt keys are present; recList/geoList
    point at existing files; filter_percent is a
    number in (0, 100]; the docking keys parse (via `load_vinalc_options`);
    recList/geoList pair up correctly (via `load_docking_targets`); at least
    one tranche was loaded; every tranche's curl_script exists and is
    unique; no two tranches share a label; and workdir exists.

    `ligList` is optional and unused: VinaLC's per-tranche ligand list is
    generated during the unpack stage (`pipeline.unpack_tranche`) from each
    tranche's downloaded ligands.

    Args:
        setup: Parsed setup dict, as returned by `load_setup`.
        tranches: Parsed tranche list, as returned by `load_tranches_tsv`.
        workdir: The run's output directory; must already exist.

    Raises:
        ValueError: For any structurally invalid input (see checks above).
        FileNotFoundError: For any referenced path that doesn't exist.
    """
    required = ["recList", "geoList", "filter_percent", "cores"]
    missing = [k for k in required if k not in setup]
    if missing:
        raise ValueError(f"setup missing required keys: {missing}")

    # Validate paths (ligList is unused, so a stale entry for it is harmless)
    for name in ("recList", "geoList"):
        p = Path(setup[name]).expanduser().resolve()
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

    load_vinalc_options(setup)

    # Fail fast on malformed recList/geoList pairing rather than during docking.
    load_docking_targets(setup)

    if not tranches:
        raise ValueError("No tranches loaded from mapping file")

    # Validate tranche scripts + uniqueness
    seen_labels: set[str] = set()
    seen_scripts: set[Path] = set()
    for t in tranches:
        if t.curl_script in seen_scripts:
            raise ValueError(f"Duplicate curl_script in mapping: {t.curl_script}")
        seen_scripts.add(t.curl_script)

        if not t.curl_script.exists():
            raise FileNotFoundError(f"curl script not found: {t.curl_script}")

        label = tranche_label(t)
        if label in seen_labels:
            raise ValueError(f"Duplicate tranche in mapping: {label}")
        seen_labels.add(label)

    # Workdir sanity
    if not workdir.exists():
        raise FileNotFoundError(f"workdir does not exist: {workdir}")

def format_summary(setup: Dict[str, str], tranches: List[Tranche], workdir: Path) -> str:
    """Render a short human-readable summary of a validated run configuration.

    Shown to the user right before Stage 0 starts (CLI) or as the result of
    clicking "Validate" (GUI), so they can sanity-check the run before any
    files are downloaded.
    """
    opts = load_vinalc_options(setup)
    logps = sorted({t.log_p for t in tranches})
    sizes = sorted({tranche_size(t) for t in tranches})
    return (
        "Inputs look valid.\n"
        f"Workdir: {workdir}\n"
        f"Tranches: {len(tranches)}\n"
        f"log_p bins ({len(logps)}): {logps}\n"
        f"Size bins ({len(sizes)}): {sizes}\n"
        f"filter_percent: {setup['filter_percent']}\n"
        f"cores: {setup['cores']}, memory_gb: {setup.get('memory_gb', 'unlimited')} -> mpirun -np {opts.mpi_ranks} "
        f"(1 master + {opts.mpi_ranks - 1} workers x ~{opts.exhaustiveness} threads)\n"
        f"vinalc options: {' '.join(opts.cli_args())}\n"
        f"recList: {setup['recList']}\n"
        f"geoList: {setup['geoList']}\n"
    )
