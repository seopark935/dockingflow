#!/usr/bin/env python3
# Usage: python3 pipeline.py --map tranches.tsv --setup setup.txt --workdir runs/run1

import argparse
from pathlib import Path

from io_parse import load_setup, load_tranches_tsv, validate_inputs, format_summary

def main() -> int:
    parser = argparse.ArgumentParser(description="Docking pipeline orchestrator")
    parser.add_argument("--map", required=True, help="Path to tranches.tsv")
    parser.add_argument("--setup", required=True, help="Path to setup.txt (KEY=VALUE)")
    parser.add_argument("--workdir", required=True, help="Output directory for this run")
    args = parser.parse_args()

    workdir = Path(args.workdir).expanduser().resolve()
    workdir.mkdir(parents=True, exist_ok=True)

    setup = load_setup(Path(args.setup))
    tranches = load_tranches_tsv(Path(args.map))

    validate_inputs(setup, tranches, workdir)

    print(format_summary(setup, tranches, workdir))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
