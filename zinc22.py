"""Pick ZINC22 3D tranches and build the pipeline's tranche map, without the CartBlanche website.

CartBlanche22's own "Download" button has problems for docking:

  1. It silently produces an *empty* script for tranches that have no 3D
     files — about 40% of what its 3D browser shows, including everything
     with 30+ heavy atoms.
  2. Its curl commands give every archive of a tranche the same output name
     (`-o H05/H05M000O.pdbqt.tgz` for generations a, d, g, q, ...), so each
     download overwrites the previous one: ~99% of the data would be lost.
  3. Its files.docking.org links often fail: generation "n" (about half the
     files) needs a login (HTTP 401), and some files are missing (404).

ZINC also publishes the same files in a public AWS S3 bucket (`zinc3d`),
readable anonymously over plain HTTPS: no AWS account or CLI needed. In a
random sample it served 92% of listed files vs. 60% from files.docking.org,
including generation "n" without a login. So this module asks CartBlanche22's
API only for the *list* of archives, and writes one curl command per archive
that downloads it from S3 (falling back to files.docking.org) to its own
unique path. The result is the same per-tranche scripts +
`curl_script tranche` map that `zinc_split.py` writes.
"""
from __future__ import annotations

import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional

import io_parse
import zinc_split

S3_BASE = "https://zinc3d.s3.amazonaws.com/"
DOCKING_ORG_BASE = "https://files.docking.org/zinc22/"
INDEX_URL = "https://cartblanche22.docking.org/tranches/get3d"
DOWNLOAD_URL = "https://cartblanche22.docking.org/tranches/3d/download"
INDEX_MAX_AGE_S = 7 * 24 * 3600
BATCH = 5000  # tranches per download-list request

# 3D files only exist for these heavy atom counts, even though CartBlanche
# lists tranches up to H49.
HAC_MIN, HAC_MAX = 4, 29

CHARGES = {"J": -4, "K": -3, "L": -2, "M": -1, "N": 0, "O": 1, "P": 2, "Q": 3, "R": 4}

# logP bins ZINC22 uses (the tranche code's P/M number / 100).
LOGP_BINS = [-5.0, -4.0, -3.0, -2.0, -1.0] + [round(x / 10, 1) for x in range(0, 51)] + [6.0, 7.0, 8.0, 9.0]

# Common screening libraries, expressed in ZINC22's axes. ZINC22 bins by heavy
# (non-hydrogen) atom count rather than molecular weight; drug-like molecules
# average ~13-14 Da per heavy atom, so e.g. 25 heavy atoms ~ 350 Da. 3D files
# stop at 29 heavy atoms (~400 Da), so wider definitions are capped there.
LIGAND_PRESETS = {
    "Quick test (~1,000 molecules)": dict(
        hac=(10, 10), logp=(3.1, 3.1), charges="Neutral only",
        help="One small tranche of about 1,000 molecules: a fast end-to-end test of your setup."),
    "Fragments (Rule of Three)": dict(
        hac=(8, 19), logp=(-1.0, 3.0), charges="Neutral only",
        help="Small, simple molecules (up to ~250 Da, logP <= 3) for fragment-based screening: "
             "weak binders that are good starting points to grow from."),
    "Lead-like": dict(
        hac=(17, 25), logp=(-1.0, 3.5), charges="Neutral and +/-1",
        help="~250-350 Da, logP <= 3.5: room to add potency and still stay drug-like. "
             "The usual choice for a first large virtual screen."),
    "Drug-like (Lipinski)": dict(
        hac=(17, 29), logp=(-1.0, 5.0), charges="Neutral and +/-1",
        help="Lipinski's rule of five (logP <= 5; MW <= 500, capped at ~400 Da here because ZINC22's 3D "
             "files stop at 29 heavy atoms). Broadest drug-like set; very large."),
}
CHARGE_PRESETS = {
    "Neutral only": ["N"],
    "Neutral and +/-1": ["M", "N", "O"],
    "All charges (-4 to +4)": list(CHARGES),
}

# An archive's key, e.g. "zinc-22a/H05/H05M000/a/H05M000-O-aaaaaa.pdbqt.tgz", from any of
# CartBlanche's formats: a files.docking.org URL (curl/wget), an s3:// URI (AWS), or an S3 URL.
_KEY_RE = re.compile(r"(?:files\.docking\.org/zinc22/|s3://zinc3d/|zinc3d\.s3\.amazonaws\.com/)(zinc-22\w/\S+)")


def _http(url: str, data: Optional[bytes] = None, headers: Optional[dict] = None, timeout: int = 300) -> bytes:
    req = urllib.request.Request(url, data=data, headers=headers or {}, method="POST" if data else "GET")
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def load_index(cache_path: Path) -> List[Dict[str, Any]]:
    """CartBlanche22's 3D tranche index, as compact dicts; cached on disk for a week.

    Each entry: {"name": "aH17P050N", "gen": "a", "hac": 17, "logp": 0.5,
    "charge": "N", "code": "H17P050", "molecules": 1229}.
    """
    if cache_path.exists() and time.time() - cache_path.stat().st_mtime < INDEX_MAX_AGE_S:
        return json.loads(cache_path.read_text())

    raw = json.loads(_http(INDEX_URL, timeout=120))
    index = []
    for t in raw["tranches"]:
        code = t["h_num"] + t["p_num"]
        hac, logp = io_parse.parse_zinc22_code(code)
        index.append({
            "name": t["generation"] + code + t["charge"],
            "gen": t["generation"], "hac": hac, "logp": logp, "charge": t["charge"],
            "code": code, "molecules": int(t.get("sum") or 0),
        })
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(index))
    return index


def select(index: List[Dict[str, Any]], hac_min: int, hac_max: int, logp_min: float, logp_max: float,
           charges: List[str]) -> List[Dict[str, Any]]:
    """Index entries within the chosen heavy-atom/logP ranges and charges."""
    return [
        t for t in index
        if hac_min <= t["hac"] <= hac_max and logp_min <= t["logp"] <= logp_max and t["charge"] in charges
    ]


def summarize(selected: List[Dict[str, Any]]) -> Dict[str, Any]:
    return {
        "tranches": len({t["code"] for t in selected}),
        "subsets": len(selected),
        "molecules": sum(t["molecules"] for t in selected),
    }


def download_command(key: str) -> str:
    """curl command fetching one archive from S3 (else files.docking.org) to `key`'s own path."""
    fetch = "curl -sS --fail --retry 3 --retry-delay 2 --create-dirs -o"
    return f"{fetch} {key} {S3_BASE}{key} || {fetch} {key} {DOCKING_ORG_BASE}{key}"


def fix_download_line(line: str) -> Optional[str]:
    """Turn one line of a CartBlanche22 downloader (curl, wget or AWS format) into `download_command`.

    Returns None for lines that don't name a ZINC22 archive.
    """
    m = _KEY_RE.search(line)
    return download_command(m.group(1)) if m else None


def build_tranche_map(selected: List[Dict[str, Any]], out_dir: Path, map_path: Path) -> Dict[str, Any]:
    """Fetch the download list for `selected`, fix it, and write per-tranche scripts + the map.

    Raises:
        ValueError: If nothing is selected, or none of it has 3D pdbqt files.
    """
    if not selected:
        raise ValueError("No tranches match these choices.")
    lines: List[str] = []
    names = [t["name"] for t in selected]
    for i in range(0, len(names), BATCH):
        body = urllib.parse.urlencode({
            "format": "pdbqt.tgz", "method": "aws", "tranches": " ".join(names[i:i + BATCH]),
        }).encode()
        resp = json.loads(_http(DOWNLOAD_URL, data=body,
                                headers={"Content-Type": "application/x-www-form-urlencoded"}))
        lines += (resp.get("data") or "").splitlines()

    fixed = [f for f in (fix_download_line(ln) for ln in lines) if f]
    if not fixed:
        raise ValueError(
            "ZINC22 has no 3D (pdbqt) files for these tranches. Try a heavy-atom range "
            f"within {HAC_MIN}-{HAC_MAX}, or more charges/logP values."
        )
    counts = zinc_split.write_tranche_scripts(fixed, "ZINC22 (via DockingFlow)", out_dir, map_path)
    return {"tranches": len(counts), "files": len(fixed)}
