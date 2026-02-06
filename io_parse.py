from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

@dataclass(frozen=True)
class Tranche:
    curl_script: Path
    log_p: float
    molecular_weight: int

def load_setup(setup_path: Path) -> Dict[str, str]:
    """Parse setup.txt with KEY=VALUE per line."""
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

def load_tranches_tsv(tsv_path: Path) -> List[Tranche]:
    """Parse tranches.tsv with header: curl_script<TAB>log_p<TAB>molecular_weight."""
    tsv_path = tsv_path.expanduser().resolve()
    if not tsv_path.exists():
        raise FileNotFoundError(f"mapping file not found: {tsv_path}")

    lines = [ln.rstrip("\n") for ln in tsv_path.read_text().splitlines() if ln.strip()]
    if not lines:
        raise ValueError("mapping file is empty")

    header = lines[0].split("\t")
    expected = ["curl_script", "log_p", "molecular_weight"]
    if header != expected:
        raise ValueError(f"Bad header. Expected {expected} but got {header}")

    tranches: List[Tranche] = []
    for i, ln in enumerate(lines[1:], start=2):
        parts = ln.split("\t")
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

        curl_script = (tsv_path.parent / curl_script_s).expanduser().resolve()
        tranches.append(Tranche(curl_script=curl_script, log_p=logp, molecular_weight=mw))

    return tranches

def validate_inputs(setup: Dict[str, str], tranches: List[Tranche], workdir: Path) -> None:
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
