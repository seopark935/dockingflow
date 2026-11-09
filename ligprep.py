"""Turn 2D SDF compounds into docking-ready 3D PDBQT, with Open Babel.

ZINC22 ships 3D PDBQT files that VinaLC can dock as they are. ChEMBL (and
most other catalogues) ship 2D SDF, which must first be given 3D
coordinates, hydrogens, partial charges and AutoDock atom types. This
module drives Open Babel (`obabel`) to do that, which is therefore required
on the docking server for ChEMBL screens but not for ZINC22 ones.

Conversion, not downloading, is the slow part of a ChEMBL screen: roughly
0.1-1 s per compound per core, so compounds are converted in batches spread
over several `obabel` processes at once.

Two things Open Babel's documentation does not specify are handled
defensively rather than assumed:

  - *How molecules are delimited* in a multi-molecule PDBQT:
    `split_pdbqt_records` accepts MODEL/ENDMDL blocks and bare
    ROOT...TORSDOF records alike.
  - *Whether compound names survive* the conversion: names are read from
    the SDF and re-attached by position, and a batch whose output count
    doesn't match its input is re-converted one compound at a time so a
    compound Open Babel drops can never shift the names of the rest.

`check_obabel` converts a two-compound sample and reports what the local
Open Babel actually did, so an install can be verified before a screen.
"""
from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterator, List, Optional, Tuple

import chembl

# Compounds per obabel process. Small enough that one failure costs little
# and progress is visible, large enough to amortise process startup.
BATCH = 100
DEFAULT_TIMEOUT_S = 600  # per batch

# A two-compound SDF (ethanol, benzene) used by `check_obabel`.
SAMPLE_SDF = """probe-ethanol
     RDKit          2D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END
$$$$
probe-benzene
     RDKit          2D

  6  6  0  0  0  0  0  0  0  0999 V2000
    0.0000    1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2990   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
M  END
$$$$
"""


@dataclass(frozen=True)
class PrepOptions:
    """How to convert SDF compounds to PDBQT (setup.txt's ligand-preparation keys).

    `gen3d` is Open Babel's 3D-generation effort: `fastest`, `fast`, `med`,
    `slow`, `slowest` (Open Babel's own default is `med`). `fast` is a good
    trade for a large screen. `minimize_steps` of 0 skips minimisation.
    """

    obabel_bin: str = "obabel"
    gen3d: str = "fast"
    forcefield: str = "MMFF94"
    minimize_steps: int = 0
    ph: Optional[str] = None
    jobs: int = 1
    timeout_s: int = DEFAULT_TIMEOUT_S

    def command(self, in_path: Path, out_path: Path) -> List[str]:
        """The obabel command line converting one SDF file to PDBQT."""
        cmd = [self.obabel_bin, str(in_path), "-O", str(out_path)]
        if self.ph:
            cmd += ["-p", self.ph]
        cmd += ["--gen3d", self.gen3d]
        if self.minimize_steps:
            # Open Babel's --minimize defaults to the Ghemical forcefield, so name ours explicitly.
            cmd += ["--minimize", "--ff", self.forcefield, "--steps", str(self.minimize_steps)]
        cmd += ["--partialcharge", "gasteiger"]
        return cmd


def available(opts: PrepOptions) -> Optional[str]:
    """The resolved path of the obabel binary, or None if it isn't installed."""
    return shutil.which(opts.obabel_bin)


def split_pdbqt_records(text: str) -> List[List[str]]:
    """Split a PDBQT file into one list of lines per molecule.

    Handles both ways Open Babel may separate molecules: MODEL/ENDMDL
    blocks, or records simply concatenated, in which case each ends with
    its single TORSDOF line.
    """
    lines = text.splitlines()
    if any(ln.startswith("MODEL") for ln in lines):
        records, current = [], None
        for ln in lines:
            if ln.startswith("MODEL"):
                current = []
            elif ln.startswith("ENDMDL"):
                if current:
                    records.append(current)
                current = None
            elif current is not None:
                current.append(ln)
        if current:
            records.append(current)
        return [r for r in records if _has_atoms(r)]

    records, current = [], []
    for ln in lines:
        current.append(ln)
        if ln.startswith("TORSDOF"):
            records.append(current)
            current = []
    if _has_atoms(current):
        records.append(current)
    return [r for r in records if _has_atoms(r)]


def _has_atoms(record: List[str]) -> bool:
    return any(ln.startswith(("ATOM", "HETATM")) for ln in record)


def _run(cmd: List[str], timeout: int, control: Any = None) -> Tuple[int, str]:
    """Run obabel, returning `(returncode, stderr)`; registers the process so a Stop can kill it."""
    proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True,
                            start_new_session=True)
    if control is not None:
        control.add_process(proc)
    try:
        _, err = proc.communicate(timeout=timeout)
        return proc.returncode, err or ""
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.communicate()
        return 1, f"timed out after {timeout}s"
    finally:
        if control is not None:
            control.discard_process(proc)


def convert_records(records: List[str], names: List[str], opts: PrepOptions,
                    control: Any = None) -> Tuple[List[Tuple[str, List[str]]], List[str]]:
    """Convert SDF records to `(name, pdbqt_lines)` pairs, keeping names aligned.

    If Open Babel returns a different number of molecules than it was
    given, the batch is redone one compound at a time, so a dropped
    compound can't shift every later name onto the wrong structure.
    """
    if not records:
        return [], []
    with tempfile.TemporaryDirectory(prefix="dockingflow_prep_") as tmp:
        tmp_dir = Path(tmp)
        sdf_path, out_path = tmp_dir / "batch.sdf", tmp_dir / "batch.pdbqt"
        sdf_path.write_text("".join(r.rstrip("\n") + "\n$$$$\n" for r in records), encoding="utf-8")
        rc, err = _run(opts.command(sdf_path, out_path), opts.timeout_s, control)
        text = out_path.read_text(errors="replace") if out_path.exists() else ""
        parsed = split_pdbqt_records(text)
        if parsed and len(parsed) == len(records):
            return list(zip(names, parsed)), []
        if len(records) == 1:
            reason = err.strip().splitlines()[-1] if err.strip() else f"obabel wrote no structure (rc={rc})"
            return [], [f"{names[0]}: {reason}"]

    # Counts disagree: fall back so every structure keeps its own name.
    converted: List[Tuple[str, List[str]]] = []
    warnings: List[str] = []
    for record, name in zip(records, names):
        if control is not None:
            control.check()
        got, warn = convert_records([record], [name], opts, control)
        converted += got
        warnings += warn
    return converted, warnings


def iter_sdf_batches(path: Path, batch: int = BATCH) -> Iterator[Tuple[List[str], List[str]]]:
    """Yield `(records, names)` batches from an SDF file (`.sdf` or `.sdf.gz`)."""
    import gzip

    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as f:  # type: ignore[operator]
        records: List[str] = []
        names: List[str] = []
        current: List[str] = []
        for line in f:
            if line.startswith("$$$$"):
                record = "".join(current)
                if record.strip():
                    records.append(record)
                    names.append(chembl.sdf_record_name(record) or f"{path.stem}_{len(records)}")
                current = []
                if len(records) >= batch:
                    yield records, names
                    records, names = [], []
            else:
                current.append(line)
        if "".join(current).strip():
            records.append("".join(current))
            names.append(chembl.sdf_record_name("".join(current)) or f"{path.stem}_{len(records)}")
        if records:
            yield records, names


def convert_sdf_file(path: Path, out_path: Path, opts: PrepOptions, control: Any = None,
                     progress: Optional[Callable[[int], None]] = None) -> Tuple[int, List[str]]:
    """Convert one SDF file to a multi-molecule PDBQT at `out_path`.

    Each molecule is written as a MODEL block carrying a `REMARK  Name =`
    line, the same shape the unpack stage expects of ZINC22 input.

    Returns `(molecules_written, warnings)`.
    """
    batches = list(iter_sdf_batches(path))
    if not batches:
        return 0, [f"{path.name}: no SDF records found"]

    warnings: List[str] = []
    written = 0
    with open(out_path, "w", encoding="utf-8") as out:
        with ThreadPoolExecutor(max_workers=max(1, opts.jobs)) as pool:
            for converted, warns in pool.map(
                    lambda b: convert_records(b[0], b[1], opts, control), batches):
                if control is not None:
                    control.check()
                warnings += warns
                for name, lines in converted:
                    written += 1
                    out.write(f"MODEL {written}\n{NAME_REMARK}{name}\n")
                    out.write("\n".join(lines) + "\nENDMDL\n")
                if progress:
                    progress(written)
    return written, warnings


NAME_REMARK = "REMARK  Name = "


def check_obabel(opts: PrepOptions) -> dict:
    """Convert a two-compound sample and report what the local Open Babel did.

    Open Babel's docs don't state how it delimits molecules in a
    multi-molecule PDBQT, nor whether compound names survive, so this
    reports both from an actual run — worth doing once on a new server
    before a ChEMBL screen.
    """
    path = available(opts)
    if not path:
        return {"ok": False, "message": (
            f"Open Babel ({opts.obabel_bin}) is not installed or not on PATH. ChEMBL screens need it to "
            "turn 2D compounds into dockable 3D structures; ZINC22 screens don't. Install it "
            "(e.g. 'conda install -c conda-forge openbabel', or ask your admin for the 'openbabel' package), "
            "or give its full path in the Ligands tab.")}
    with tempfile.TemporaryDirectory(prefix="dockingflow_check_") as tmp:
        sdf = Path(tmp) / "sample.sdf"
        sdf.write_text(SAMPLE_SDF)
        records = list(chembl.split_sdf_records(SAMPLE_SDF))
        names = [chembl.sdf_record_name(r) or "?" for r in records]
        out = Path(tmp) / "sample.pdbqt"
        rc, err = _run(opts.command(sdf, out), min(opts.timeout_s, 120))
        text = out.read_text(errors="replace") if out.exists() else ""
        parsed = split_pdbqt_records(text)
        if not parsed:
            tail = err.strip().splitlines()[-3:] if err.strip() else []
            return {"ok": False, "message": (
                f"Open Babel at {path} ran but produced no structures (rc={rc}). "
                + (" / ".join(tail) if tail else "It printed no error.")
                + f"\nCommand: {' '.join(opts.command(Path('in.sdf'), Path('out.pdbqt')))}")}
        has_3d = any(_looks_3d(r) for r in parsed)
        return {"ok": True, "version": _version(opts), "path": path, "molecules": len(parsed),
                "expected": len(records), "delimiter": "MODEL/ENDMDL" if "MODEL" in text else "concatenated",
                "names_preserved": all(n in text for n in names), "has_3d": has_3d,
                "message": (
                    f"Open Babel {_version(opts)} at {path} converted {len(parsed)} of {len(records)} sample "
                    f"compounds; 3D coordinates {'present' if has_3d else 'MISSING'}. Ready for ChEMBL screens."
                    if len(parsed) == len(records) and has_3d else
                    f"Open Babel at {path} converted {len(parsed)} of {len(records)} sample compounds and 3D "
                    f"coordinates are {'present' if has_3d else 'missing'}. Check the --gen3d setting.")}


def _looks_3d(record: List[str]) -> bool:
    """Whether a PDBQT record's atoms have a non-zero z coordinate (i.e. real 3D, not flattened 2D)."""
    for line in record:
        if line.startswith(("ATOM", "HETATM")):
            try:
                if abs(float(line[46:54])) > 1e-6:
                    return True
            except ValueError:
                continue
    return False


def _version(opts: PrepOptions) -> str:
    try:
        out = subprocess.run([opts.obabel_bin, "-V"], capture_output=True, text=True, timeout=30)
        return (out.stdout or out.stderr).strip().splitlines()[0] if (out.stdout or out.stderr) else "?"
    except (OSError, subprocess.SubprocessError):
        return "?"
