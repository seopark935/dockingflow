#!/usr/bin/env python3
"""Mock VinaLC binary for offline testing.

Accepts VinaLC's real command line (`--recList --ligList --geoList` plus
tuning flags, which are ignored) and writes output in VinaLC's real format
and location: one gzipped `<recList>_<ligList>.pdbqt.gz` (plus `.log.gz`)
in the working directory, containing a record per (receptor, ligand) job::

    REMARK RECEPTOR <receptor path>
    REMARK LIGAND LIGAND <n>
    MODEL 1
    REMARK VINA RESULT: ...
    <ligand lines, including its REMARK Name>
    ENDMDL
    ...

Like the real binary, ligands are only the MODEL...ENDMDL blocks of each
ligand file, numbered per receptor across all ligand files in order, and
records are written out of input order (VinaLC writes in completion order),
so the pipeline's ranking has to map `LIGAND <n>` back via its own index.
Affinities are deterministic fakes derived from the ligand text.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
from pathlib import Path


def fake_affinity(key: str) -> float:
    digest = hashlib.sha256(key.encode()).hexdigest()
    return -(4.0 + (int(digest, 16) % 6000) / 1000.0)  # roughly -4.0 to -10.0


def read_models(path: Path) -> list[str]:
    """MODEL...ENDMDL block bodies, parsed the way VinaLC's master rank does."""
    models, current = [], None
    for line in path.read_text().splitlines():
        if line.startswith("MODEL"):
            current = []
        elif line.startswith("ENDMDL"):
            if current is not None:
                models.append("\n".join(current))
            current = None
        elif current is not None:
            current.append(line)
    return models


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--recList", required=True)
    parser.add_argument("--ligList", required=True)
    parser.add_argument("--geoList", required=True)
    args, _unknown = parser.parse_known_args()

    receptors = [ln.strip() for ln in Path(args.recList).read_text().splitlines() if ln.strip()]
    geos = [ln for ln in Path(args.geoList).read_text().splitlines() if ln.strip()]
    if len(receptors) != len(geos):
        print("Error: Receptor and geometry lists are not equal")
        return 1
    lig_files = [Path(p.strip()) for p in Path(args.ligList).read_text().splitlines() if p.strip()]

    records = []
    for receptor in receptors:
        n = 0
        for lig_file in lig_files:
            for body in read_models(lig_file):
                n += 1
                best = fake_affinity(receptor + body)
                poses = "".join(
                    f"MODEL {i}\nREMARK VINA RESULT: {best + (i - 1) * 0.3:9.1f}      0.000      0.000\n{body}\nENDMDL\n"
                    for i in range(1, 4)
                )
                records.append(f"REMARK RECEPTOR {receptor}\nREMARK LIGAND LIGAND {n}\n{poses}\n")

    stem = f"{args.recList}_{args.ligList}"
    with gzip.open(f"{stem}.pdbqt.gz", "wt") as out:
        out.write("".join(reversed(records)))  # not input order, like VinaLC
    with gzip.open(f"{stem}.log.gz", "wt") as log:
        log.write(f"mock vinalc: {len(records)} job(s)\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
