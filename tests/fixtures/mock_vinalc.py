#!/usr/bin/env python3
"""Mock vinalc binary for offline testing.

Reads a vinalc-style config file (KEY = VALUE lines) and, for every ligand in
its ligand_list, writes a fake docked *_out.pdbqt file with deterministic
'REMARK VINA RESULT:' lines so downstream parsing can be exercised without a
real docking binary.
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


def parse_config(path: Path) -> dict[str, str]:
    cfg: dict[str, str] = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or "=" not in line:
            continue
        key, val = line.split("=", 1)
        cfg[key.strip()] = val.strip()
    return cfg


def fake_affinity(ligand_name: str) -> float:
    digest = hashlib.sha256(ligand_name.encode()).hexdigest()
    return -(4.0 + (int(digest, 16) % 6000) / 1000.0)  # roughly -4.0 to -10.0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    cfg = parse_config(Path(args.config))
    ligand_list = Path(cfg["ligand_list"])
    out_dir = Path(cfg["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    ligands = [Path(p) for p in ligand_list.read_text().splitlines() if p.strip()]
    for ligand in ligands:
        best = fake_affinity(ligand.stem)
        modes = "\n".join(
            f"MODEL {i}\nREMARK VINA RESULT:    {best + (i - 1) * 0.3:7.2f}      0.000      0.000\nENDMDL"
            for i in range(1, 4)
        )
        (out_dir / f"{ligand.stem}_out.pdbqt").write_text(
            f"REMARK  Fake docking output for {ligand.name}\n{modes}\n"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
