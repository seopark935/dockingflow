"""Pick ChEMBL compounds and build the pipeline's tranche map, as an alternative to ZINC22.

ZINC22 ships ready-to-dock 3D PDBQT files; ChEMBL does not. ChEMBL serves
2D SDF, so every compound needs 3D coordinates and AutoDock atom typing
before docking. The unpack stage does that with Open Babel
(`pipeline.sdf_to_pdbqt`), which must be installed on the docking server —
see `server_check.sh`. That conversion, not the download, is the slow part
of a ChEMBL screen.

Two ways to fetch compounds:

  - `source="api"` (default): ChEMBL's REST API, filtered by molecular
    weight and logP, requested as SDF in pages of 1000. The selection is
    split into molecular-weight bins, one tranche each, so a screen is
    resumable and results appear bin by bin. Allow roughly 10 s per 1000
    compounds.
  - `source="ftp"`: the whole release as one ~1 GB file from EBI's FTP
    site. Far quicker to fetch than paging ~2 million compounds through
    the API, but it arrives as a single tranche: no per-tranche
    checkpoints, and one very long docking run.

Compound identity: a ChEMBL SDF record has an empty title line and carries
its identifier in a `> <chembl_id>` data field, so `sdf_record_name` reads
that rather than the title.
"""
from __future__ import annotations

import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Tuple

import zinc_split

API_URL = "https://www.ebi.ac.uk/chembl/api/data/molecule"
FTP_DIR = "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/"
PAGE_SIZE = 1000  # the API's maximum
COUNT_TIMEOUT_S = 180  # a filtered count can take a couple of minutes

# Molecular weight range the bins span when a preset leaves one end open.
MW_FLOOR, MW_CEILING = 0.0, 1000.0
DEFAULT_BIN_WIDTH = 25.0  # Da per tranche

# Screening libraries in ChEMBL's axes (molecular weight in Da, Crippen
# logP as ChEMBL's "alogp"). Counts in the help text are from ChEMBL 37 and
# shift between releases, so the GUI always shows a live count.
LIGAND_PRESETS = {
    "Quick test (~900 molecules)": dict(
        mw=(200.0, 201.0), logp=(None, None),
        help="One narrow weight window, about 900 compounds: a fast end-to-end test of your setup."),
    "Fragments (Rule of Three)": dict(
        mw=(None, 300.0), logp=(None, 3.0),
        help="Small, simple compounds (<= 300 Da, logP <= 3) for fragment-based screening: weak binders "
             "that are good starting points to grow from. About 230,000 compounds."),
    "Lead-like": dict(
        mw=(250.0, 350.0), logp=(None, 3.5),
        help="250-350 Da, logP <= 3.5: room to add potency and still stay drug-like. The usual choice "
             "for a first large virtual screen. About 370,000 compounds."),
    "Drug-like (Lipinski)": dict(
        mw=(None, 500.0), logp=(None, 5.0),
        help="Lipinski's rule of five (<= 500 Da, logP <= 5). The broadest drug-like set, over a million "
             "compounds; expect a long screen."),
    "Whole library (every small molecule)": dict(
        mw=(None, None), logp=(None, None),
        help="Every small molecule in ChEMBL, about 1.9 million compounds. This is the set used for "
             "large published screens; expect weeks of docking."),
}

_ID_TAG_RE = re.compile(r"^>\s*<chembl_id>", re.M)


# ChEMBL's API returns an occasional HTTP 500 under load, so every read is
# retried; the generated download scripts pass curl --retry for the same reason.
RETRIES, RETRY_WAIT_S = 4, 3


def _http(url: str, timeout: int = 120) -> bytes:
    last: Exception | None = None
    for attempt in range(RETRIES):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "DockingFlow"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, OSError) as exc:
            status = getattr(exc, "code", None)
            if status is not None and status < 500 and status != 429:
                raise  # a real rejection (bad filter, 404): retrying won't help
            last = exc
            if attempt < RETRIES - 1:
                time.sleep(RETRY_WAIT_S * (attempt + 1))
    raise RuntimeError(f"ChEMBL did not answer after {RETRIES} attempts ({last}). Try again in a few minutes.")


def _filter_params(mw_min: Optional[float], mw_max: Optional[float],
                   logp_min: Optional[float], logp_max: Optional[float]) -> List[Tuple[str, str]]:
    """ChEMBL API query parameters for a molecular weight / logP window.

    Entries without a structure are excluded. They can't be docked, and
    worse, a single one of them makes the SDF endpoint fail its *whole*
    page with "Molecule has no structure records" (HTTP 404), so leaving
    them in would lose up to 1000 good compounds at a time. Counting uses
    the same filter, so the count matches what actually downloads.
    """
    params: List[Tuple[str, str]] = [
        ("molecule_type", "Small molecule"),
        ("molecule_structures__isnull", "false"),
    ]
    for key, value in (("full_mwt__gte", mw_min), ("full_mwt__lte", mw_max),
                       ("alogp__gte", logp_min), ("alogp__lte", logp_max)):
        if value is not None:
            params.append((f"molecule_properties__{key}", f"{value:g}"))
    return params


def count(mw_min: Optional[float] = None, mw_max: Optional[float] = None,
          logp_min: Optional[float] = None, logp_max: Optional[float] = None) -> int:
    """How many ChEMBL small molecules match this window (one cheap API call)."""
    params = _filter_params(mw_min, mw_max, logp_min, logp_max) + [("limit", "1"), ("format", "json")]
    url = f"{API_URL}?{urllib.parse.urlencode(params)}"
    return int(json.loads(_http(url, COUNT_TIMEOUT_S))["page_meta"]["total_count"])


def page_url(mw_min: Optional[float], mw_max: Optional[float], logp_min: Optional[float],
             logp_max: Optional[float], offset: int, limit: int = PAGE_SIZE) -> str:
    """URL of one page of matching compounds as SDF.

    Ordered by ChEMBL id so paging is stable: without an explicit order the
    API may return rows in a different order per request, which would
    duplicate some compounds and miss others across pages.
    """
    params = _filter_params(mw_min, mw_max, logp_min, logp_max) + [
        ("order_by", "molecule_chembl_id"), ("limit", str(limit)), ("offset", str(offset))]
    return f"{API_URL}.sdf?{urllib.parse.urlencode(params)}"


def bin_edges(mw_min: Optional[float], mw_max: Optional[float], width: float = DEFAULT_BIN_WIDTH) -> List[Tuple[float, float]]:
    """Split a weight range into `width`-Da bins, each of which becomes one tranche."""
    lo = MW_FLOOR if mw_min is None else float(mw_min)
    hi = MW_CEILING if mw_max is None else float(mw_max)
    if hi <= lo:
        raise ValueError(f"The maximum weight ({hi:g}) must be above the minimum ({lo:g}).")
    edges = []
    start = lo
    while start < hi:
        end = min(start + width, hi)
        edges.append((start, end))
        start = end
    return edges


def bin_label(mw_lo: float, mw_hi: float) -> str:
    """Tranche label for a weight bin, e.g. `MW250-275` (a directory name, so no spaces or symbols)."""
    return f"MW{mw_lo:g}-{mw_hi:g}".replace(".", "p")


def summarize(mw_min: Optional[float], mw_max: Optional[float], logp_min: Optional[float],
              logp_max: Optional[float], source: str = "api") -> Dict[str, Any]:
    """Compound and tranche counts for the current selection, for the GUI's live summary."""
    if source == "ftp":
        return {"tranches": 1, "molecules": count(), "whole_release": True}
    n = count(mw_min, mw_max, logp_min, logp_max)
    return {"tranches": len(bin_edges(mw_min, mw_max)), "molecules": n, "whole_release": False}


def release_file() -> str:
    """Name of the current release's full SDF on EBI's FTP site (e.g. `chembl_37.sdf.gz`)."""
    listing = _http(FTP_DIR, timeout=60).decode("utf-8", errors="replace")
    names = sorted(set(re.findall(r'href="(chembl_\d+\.sdf\.gz)"', listing)))
    if not names:
        raise ValueError(f"No chembl_*.sdf.gz found at {FTP_DIR}")
    return names[-1]


def build_tranche_map(out_dir: Path, map_path: Path, mw_min: Optional[float] = None,
                      mw_max: Optional[float] = None, logp_min: Optional[float] = None,
                      logp_max: Optional[float] = None, source: str = "api",
                      width: float = DEFAULT_BIN_WIDTH) -> Dict[str, Any]:
    """Write one download script per tranche plus the `curl_script tranche` map.

    With `source="api"` each weight bin becomes a tranche whose script
    fetches its compounds as paged SDF; empty bins are skipped. With
    `source="ftp"` there is a single tranche that fetches the whole
    release.

    Raises:
        ValueError: If nothing matches the selection.
    """
    if source == "ftp":
        name = release_file()
        line = (f"curl -sS --fail --retry 3 --retry-delay 2 --create-dirs -o {name} {FTP_DIR}{name}")
        label = name.split(".")[0].upper().replace("_", "")  # chembl_37.sdf.gz -> CHEMBL37
        counts = _write_scripts({label: [line]}, f"ChEMBL release {name}", out_dir, map_path)
        return {"tranches": 1, "files": 1, "molecules": None, "labels": list(counts)}

    groups: Dict[str, List[str]] = {}
    total = 0
    for lo, hi in bin_edges(mw_min, mw_max, width):
        # Bins meet end to end, so the upper edge is exclusive except on the
        # last bin; otherwise a compound exactly on an edge lands in both.
        n = count(lo, hi, logp_min, logp_max)
        if not n:
            continue
        total += n
        label = bin_label(lo, hi)
        groups[label] = [
            "curl -sS --fail --retry 3 --retry-delay 2 --create-dirs "
            f"-o page_{page + 1:05d}.sdf '{page_url(lo, hi, logp_min, logp_max, page * PAGE_SIZE)}'"
            for page in range((n + PAGE_SIZE - 1) // PAGE_SIZE)
        ]
    if not groups:
        raise ValueError("No ChEMBL compounds match these filters. Widen the weight or logP range.")

    counts = _write_scripts(groups, "ChEMBL API (via DockingFlow)", out_dir, map_path)
    return {"tranches": len(counts), "files": sum(counts.values()), "molecules": total, "labels": list(counts)}


def _write_scripts(groups: Dict[str, List[str]], source: str, out_dir: Path, map_path: Path) -> Dict[str, int]:
    """Write a download script per label plus the map file (the ChEMBL twin of `zinc_split.write_tranche_scripts`)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    for old in out_dir.glob("*.curl"):
        old.unlink()
    rows = []
    import os

    for label, cmds in groups.items():
        script = out_dir / f"{label}.curl"
        script.write_text(zinc_split.SCRIPT_HEADER.format(source=source, code=label)
                          + "\n".join(cmds) + "\n" + zinc_split.SCRIPT_FOOTER)
        rows.append(f"{os.path.relpath(script, map_path.parent)}\t{label}")
    map_path.parent.mkdir(parents=True, exist_ok=True)
    map_path.write_text("curl_script\ttranche\n" + "\n".join(rows) + "\n")
    return {label: len(cmds) for label, cmds in groups.items()}


def split_sdf_records(text: str) -> Iterator[str]:
    """Yield each record of an SDF file, without its `$$$$` terminator."""
    for record in text.split("$$$$"):
        if record.strip():
            # Exactly one newline: the one after "$$$$". Stripping more would
            # eat an SDF's (often empty) title line and shift every line up.
            yield record[1:] if record.startswith("\n") else record


def sdf_record_name(record: str) -> Optional[str]:
    """A record's compound id: its `> <chembl_id>` field, else its title line.

    ChEMBL's SDF leaves the title line empty and puts the identifier in a
    data field, while SDFs from other sources usually use the title.
    """
    lines = record.splitlines()
    for i, line in enumerate(lines):
        if _ID_TAG_RE.match(line) and i + 1 < len(lines) and lines[i + 1].strip():
            return lines[i + 1].strip()
    return lines[0].strip() if lines and lines[0].strip() else None
