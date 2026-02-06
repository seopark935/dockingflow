#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
from datetime import datetime
from pathlib import Path

from io_parse import load_setup, load_tranches_tsv, validate_inputs, format_summary, Tranche


def tranche_label(t: Tranche) -> str:
    """Filesystem-safe, deterministic tranche label."""
    # Normalize logp for folder name: 5.0 -> "5", 0.50 -> "0.5"
    logp = f"{t.log_p:.2f}".rstrip("0").rstrip(".")
    return f"LP{logp}_MW{t.molecular_weight}"


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def materialize_tranche(workdir: Path, tranche: Tranche, setup_path: Path) -> Path:
    """
    Create tranche workspace folder and snapshot key inputs.
    Resumable: does not overwrite existing status/meta.
    """
    tdir = workdir / "tranches" / tranche_label(tranche)
    inputs_dir = tdir / "inputs"
    inputs_dir.mkdir(parents=True, exist_ok=True)

    # Snapshot artifacts for audit trail
    shutil.copy2(tranche.curl_script, inputs_dir / "curl_script.curl")
    shutil.copy2(setup_path, inputs_dir / "setup_used.txt")

    # Status file: created once
    status_path = tdir / "status.txt"
    if not status_path.exists():
        write_text(status_path, "INIT\n")

    # Meta file: created once
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


def main() -> int:
    parser = argparse.ArgumentParser(description="Docking pipeline orchestrator (stage 0)")
    parser.add_argument("--map", required=True, help="Path to tranches.tsv (tab-separated)")
    parser.add_argument("--setup", required=True, help="Path to setup.txt (KEY=VALUE)")
    parser.add_argument("--workdir", required=True, help="Output directory for this run")
    args = parser.parse_args()

    workdir = Path(args.workdir).expanduser().resolve()
    workdir.mkdir(parents=True, exist_ok=True)

    setup_path = Path(args.setup).expanduser().resolve()
    map_path = Path(args.map).expanduser().resolve()

    setup = load_setup(setup_path)
    tranches = load_tranches_tsv(map_path)

    validate_inputs(setup, tranches, workdir)

    print(format_summary(setup, tranches, workdir))

    (workdir / "tranches").mkdir(exist_ok=True)

    created = 0
    for t in tranches:
        materialize_tranche(workdir, t, setup_path)
        created += 1

    print(f"\nStage 0 complete: prepared {created} tranche workspaces in:")
    print(f"  {workdir / 'tranches'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
