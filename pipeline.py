#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
import subprocess
import time
from datetime import datetime
from pathlib import Path

from io_parse import (
    Tranche,
    format_summary,
    load_setup,
    load_tranches_tsv,
    validate_inputs,
)

# ----------------------------
# Small utilities
# ----------------------------

FINAL_STATUSES = {"DOWNLOADED", "DOCKED", "DONE"}


def tranche_label(t: Tranche) -> str:
    """Filesystem-safe, deterministic tranche label."""
    logp = f"{t.log_p:.2f}".rstrip("0").rstrip(".")
    return f"LP{logp}_MW{t.molecular_weight}"


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def read_status(tdir: Path) -> str:
    p = tdir / "status.txt"
    if not p.exists():
        return "INIT"
    return p.read_text(encoding="utf-8").strip() or "INIT"


def write_status(tdir: Path, status: str) -> None:
    write_text(tdir / "status.txt", status + "\n")


def run_cmd(cmd: list[str], cwd: Path, stdout_path: Path, stderr_path: Path) -> tuple[int, float]:
    """Run a command, redirecting stdout/stderr to files. Returns (rc, seconds)."""
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    stderr_path.parent.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    with open(stdout_path, "w", encoding="utf-8") as out, open(stderr_path, "w", encoding="utf-8") as err:
        proc = subprocess.Popen(cmd, cwd=str(cwd), stdout=out, stderr=err, text=True)
        rc = proc.wait()
    return rc, time.time() - t0


def count_pdbqt_gz(download_dir: Path) -> int:
    return sum(1 for _ in download_dir.rglob("*.pdbqt.gz"))


# ----------------------------
# Stage 0: workspace creation
# ----------------------------

def materialize_tranche(workdir: Path, tranche: Tranche, setup_path: Path) -> Path:
    """
    Create tranche workspace folder and snapshot key inputs.
    Resumable: does not overwrite existing status/meta.
    """
    tdir = workdir / "tranches" / tranche_label(tranche)
    inputs_dir = tdir / "inputs"
    inputs_dir.mkdir(parents=True, exist_ok=True)

    # Snapshot artifacts for audit trail (always refresh; curl scripts may update)
    shutil.copy2(tranche.curl_script, inputs_dir / "curl_script.curl")
    shutil.copy2(setup_path, inputs_dir / "setup_used.txt")

    # Status/meta created once
    status_path = tdir / "status.txt"
    if not status_path.exists():
        write_status(tdir, "INIT")

    meta_path = tdir / "meta.txt"
    if not meta_path.exists():
        write_text(
            meta_path,
            f"created_at={datetime.now().isoformat(timespec='seconds')}\n"
            f"log_p={tranche.log_p}\n"
            f"molecular_weight={tranche.molecular_weight}\n"
            f"curl_script_original={tranche.curl_script}\n"
        )

    return tdir


# ----------------------------
# Stage 1: download
# ----------------------------

def download_tranche(tdir: Path) -> None:
    status = read_status(tdir)
    if status in FINAL_STATUSES:
        print(f"[skip] {tdir.name}: status={status}")
        return

    print(f"[download] {tdir.name}")

    download_dir = tdir / "download"
    logs_dir = tdir / "logs"
    download_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    curl_script = tdir / "inputs" / "curl_script.curl"
    if not curl_script.exists():
        write_status(tdir, "FAILED_DOWNLOAD")
        raise FileNotFoundError(f"Missing curl script snapshot: {curl_script}")

    stdout_log = logs_dir / "download.stdout.log"
    stderr_log = logs_dir / "download.stderr.log"

    # Run curl script in download_dir so --create-dirs writes underneath it.
    rc, dt = run_cmd(["bash", str(curl_script)], cwd=download_dir, stdout_path=stdout_log, stderr_path=stderr_log)

    n_files = count_pdbqt_gz(download_dir)
    write_text(tdir / "download_timing.txt", f"rc={rc}\nseconds={dt:.2f}\nfiles={n_files}\n")

    if rc != 0:
        write_status(tdir, "FAILED_DOWNLOAD")
        raise RuntimeError(f"Download failed for {tdir.name} (rc={rc}). See {stderr_log}")

    if n_files == 0:
        write_status(tdir, "FAILED_DOWNLOAD")
        raise RuntimeError(f"No *.pdbqt.gz files downloaded for {tdir.name}. See logs in {logs_dir}")

    write_status(tdir, "DOWNLOADED")


# ----------------------------
# Cleaning
# ----------------------------

def clean_workdir(workdir: Path) -> None:
    """
    Delete all generated tranche directories under workdir/tranches/.
    Intended for test resets. Keeps the workdir itself.
    """
    tranches_dir = workdir / "tranches"
    if not tranches_dir.exists():
        print("Nothing to clean (no tranches directory).")
        return

    # Safety guard: never delete root or home by accident
    if tranches_dir == Path("/") or tranches_dir == Path.home():
        raise RuntimeError(f"Refusing to delete dangerous path: {tranches_dir}")

    shutil.rmtree(tranches_dir)
    print(f"Deleted generated tranche directories: {tranches_dir}")

# ----------------------------
# Main
# ----------------------------

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Docking pipeline orchestrator (Stage 0 + Download stage)"
    )
    parser.add_argument("--map", required=True, help="Path to tranches.tsv (tab-separated)")
    parser.add_argument("--setup", required=True, help="Path to setup.txt (KEY=VALUE)")
    parser.add_argument("--workdir", required=True, help="Output directory for this run")
    parser.add_argument("--clean", action="store_true", help="Delete generated outputs in workdir and reset statuses")
    parser.add_argument("--nuke", action="store_true", help="Delete the entire workdir (DANGEROUS)")
    parser.add_argument("--only-stage0", action="store_true", help="Only create tranche workspaces (no downloads)")
    args = parser.parse_args()

    workdir = Path(args.workdir).expanduser().resolve()

    # Clean/nuke modes should not require valid inputs.
    if args.nuke:
        nuke_workdir(workdir)
        return 0

    workdir.mkdir(parents=True, exist_ok=True)

    if args.clean:
        clean_workdir(workdir)
        return 0

    setup_path = Path(args.setup).expanduser().resolve()
    map_path = Path(args.map).expanduser().resolve()

    setup = load_setup(setup_path)
    tranches = load_tranches_tsv(map_path)

    validate_inputs(setup, tranches, workdir)

    print(format_summary(setup, tranches, workdir))

    (workdir / "tranches").mkdir(exist_ok=True)

    # Stage 0: create tranche directories + snapshot inputs
    tranche_dirs: list[Path] = []
    for t in tranches:
        tdir = materialize_tranche(workdir, t, setup_path)
        tranche_dirs.append(tdir)

    print(f"\nStage 0 complete: prepared {len(tranche_dirs)} tranche workspaces.")

    if args.only_stage0:
        return 0

    # Stage 1: download each tranche sequentially (safe default)
    for tdir in tranche_dirs:
        download_tranche(tdir)

    print("\nDownload stage complete for all tranches.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
