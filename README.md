# dockingflow

A resumable virtual-screening pipeline: downloads ZINC ligand tranches,
docks them against a receptor with a parallel AutoDock Vina backend
(VinaLC-style), and ranks/filters the results down to top hits. Comes with
both a CLI (`pipeline.py`) and a desktop GUI (`gui.py`) built on the same
underlying code.

- [Glossary](#glossary)
- [Quick start](#quick-start)
- [Architecture](#architecture)
- [Pipeline stages](#pipeline-stages-in-detail)
- [CLI reference](#cli-reference)
- [GUI](#gui)
- [Configuration file reference](#configuration-file-reference)
- [Directory layout](#directory-layout-produced-by-a-run)
- [Placeholders used for offline testing](#placeholders-used-for-offline-testing)
- [Testing](#testing)
- [Troubleshooting](#troubleshooting)

## Glossary

If you're not coming from a computational chemistry background, these terms
show up throughout the config files and code:

- **Docking**: computationally predicting how (and how strongly) a small
  molecule (a *ligand*) binds to a target protein (the *receptor*), by
  searching many possible orientations/poses and scoring each one.
- **Receptor**: the target protein structure, in `.pdbqt` format (PDB
  coordinates plus AutoDock-specific atom typing/charge info).
- **Ligand**: a small candidate molecule being screened against the
  receptor, also in `.pdbqt` format.
- **ZINC**: a public, free database of purchasable compounds
  (https://zinc20.docking.org), pre-organized into downloadable *tranches*.
- **Tranche**: a bucket of ZINC compounds sharing a property range — here,
  a `(log_p, molecular_weight)` bin. Screening campaigns typically process
  many tranches to cover a broad chemical space.
- **LogP**: a measure of a molecule's lipophilicity (fat solubility).
  Screens are often bucketed by LogP because it correlates with drug-like
  properties.
- **Grid box / search box**: the 3D region (a center point + x/y/z
  dimensions, in angstroms) the docking engine searches within — normally
  placed around a known or suspected binding pocket on the receptor.
- **Binding affinity**: the docking engine's predicted binding strength,
  in kcal/mol. More negative = predicted to bind more strongly. This
  pipeline ranks ligands by their best (most negative) pose.
- **AutoDock Vina / VinaLC**: Vina is the widely-used, free docking engine
  this pipeline's config format and output parsing are built around. VinaLC
  is an MPI-parallelized version of the same underlying algorithm, built
  for running large screens across many CPU cores at once — this is what
  the `cores` and `mpirun` options in this pipeline target.
- **Exhaustiveness / energy_range / num_modes**: Vina-family tuning knobs.
  `energy_range` controls how wide a window of binding energies gets
  reported per ligand; `num_modes` caps how many distinct poses are
  reported.

## Quick start

```bash
python3 pipeline.py --map tranches.txt --setup setup.txt --workdir ./run1
```

This uses the placeholder receptor/ligand data checked into this repo, so
it runs fully offline out of the box — see
[Placeholders used for offline testing](#placeholders-used-for-offline-testing)
before pointing it at a real screen.

## Architecture

```
                 ┌───────────────┐        ┌───────────────┐
                 │   pipeline.py │        │    gui.py     │
                 │  (CLI, main)  │        │ (PipelineAPI) │
                 └───────┬───────┘        └───────┬───────┘
                         │   both call the same    │
                         │      stage functions     │
                         ▼                          ▼
        materialize_tranche → download_tranche → unpack_tranche
                    → dock_tranche → rank_and_filter_tranche
                                    │
                                    ▼
                         combine_results (workdir-level)
                                    │
                    ┌───────────────┴───────────────┐
                    ▼                                ▼
             io_parse.py                     tests/fixtures/
      (parses + validates setup.txt,          mock_vinalc.py
       tranches.txt, recList/geoList)       (fake docking binary
                                                for offline tests)
```

`pipeline.py` and `gui.py` never duplicate pipeline logic — they're both
thin orchestration layers over the same stage functions in `pipeline.py`.
The one behavioral difference: the CLI's `main()` stops the whole run on
the first tranche that raises, while the GUI's runner isolates failures
per-tranche (records the error, moves on to the next tranche) so a partial
run still shows useful results in the window. See the docstrings on
`pipeline.main` and `gui.PipelineAPI._run` for the reasoning.

## Pipeline stages (in detail)

Each tranche progresses through statuses tracked in `<tranche>/status.txt`,
so a re-run skips whatever's already done:

```
INIT → DOWNLOADED → UNPACKED → DOCKED → DONE
```

(or `FAILED_<STAGE>` if that stage errors — e.g. `FAILED_DOWNLOAD`,
`FAILED_UNPACK`, `FAILED_DOCK`, `FAILED_RANK`). A failed tranche is
re-attempted from its failed stage on the next run, since `FAILED_*`
statuses don't count as "past" any stage.

1. **Stage 0 — workspace** (`materialize_tranche`): creates
   `workdir/tranches/<label>/` per tranche (label derived from its LogP/MW,
   e.g. `LP5_MW200`) and snapshots that tranche's curl script + the
   setup.txt used, for an audit trail. Safe to re-run — never overwrites an
   existing `status.txt` or `meta.txt`.
2. **Stage 1 — download** (`download_tranche`): runs the tranche's
   snapshotted curl script, then verifies at least one `*.pdbqt.gz` file
   was produced. Logs stdout/stderr to `<tranche>/logs/download.*.log`.
3. **Stage 2 — unpack** (`unpack_tranche`): decompresses every downloaded
   `*.pdbqt.gz` into a flat `<tranche>/ligands/` directory and writes
   `<tranche>/ligand_list.txt` (one ligand path per line) — this is the
   file the docking stage feeds to the docking binary.
4. **Stage 3 — dock** (`dock_tranche`): for every `(receptor, grid box)`
   pair defined by `recList.txt`/`geoList.txt`, writes a vinalc-style
   config (`build_vinalc_config`) and runs the docking binary — optionally
   wrapped in `mpirun -np <cores>` — against the tranche's full ligand
   list. Each target gets its own subdirectory under
   `<tranche>/docking/target<i>_<receptor>/` so multi-target runs never
   collide.
5. **Stage 4 — rank + filter** (`rank_and_filter_tranche`): parses the
   first `REMARK VINA RESULT:` line (the best pose — Vina/VinaLC always
   list poses best-first) out of every docked `*_out.pdbqt` file, ranks
   ligands ascending by affinity (most negative = best), and writes
   `<tranche>/results/ranked_all.tsv` (everything) and
   `<tranche>/results/top_hits.tsv` (top `filter_percent`%, at least 1
   ligand).
6. **Combine** (`combine_results`, run once per full pipeline invocation,
   not per-tranche): merges every tranche's `top_hits.tsv` into one
   workdir-level ranking at `workdir/results/top_hits_combined.tsv` — this
   is the final deliverable of a screening run.

## CLI reference

```
python3 pipeline.py --map <tranches.tsv> --setup <setup.txt> --workdir <dir> [options]
```

| Flag | Description |
|---|---|
| `--map PATH` | Required. Path to the tranche mapping TSV. |
| `--setup PATH` | Required. Path to `setup.txt`. |
| `--workdir PATH` | Required. Output directory for this run. |
| `--only-stage0` | Stop after creating tranche workspaces (no downloads). |
| `--only-download` | Stop after the download stage. |
| `--only-unpack` | Stop after unpacking downloaded ligands. |
| `--vinalc-bin PATH` | Docking binary to invoke (default `vinalc`). |
| `--mpirun-bin PATH` | MPI launcher to wrap the docking binary with (default `mpirun`). |
| `--no-mpirun` | Invoke the docking binary directly, skipping the MPI launcher. |
| `--clean` | Delete generated tranche outputs under the workdir; keeps the workdir itself. Does not require `--map`/`--setup` to be valid. |
| `--nuke` | Delete the entire workdir. **Destructive and irreversible.** Does not require `--map`/`--setup` to be valid. |

## GUI

A desktop app (not a browser tab) built with [pywebview](https://pywebview.flowrl.com/):
Python hosts the pipeline logic and exposes it to a bundled HTML/CSS/JS page
rendered in a native window.

### Setup

pywebview isn't in the standard library and this machine's Python is an
externally-managed Homebrew install, so it's kept in a project-local
virtual environment rather than installed system-wide:

```bash
python3 -m venv .venv
.venv/bin/pip install pywebview
```

### Running

```bash
.venv/bin/python3 gui.py
```

### What it does

- **Configuration panel** — set the setup file, tranches map, work
  directory, docking binary, and MPI launcher, with native file/folder
  pickers. Opens pre-filled with this repo's own `setup.txt`/`tranches.txt`.
- **Validate** — runs the same checks `pipeline.py` runs before Stage 0,
  without downloading or docking anything.
- **Run pipeline** — runs all stages on a background thread so the window
  stays responsive; the tranche table, log, and results panel update live
  (polled every ~700ms).
- **Tranche table** — one row per tranche, showing its current status as a
  colored badge, and any per-tranche error inline.
- **Log panel** — a running, auto-scrolling log of what the pipeline is
  doing.
- **Top hits table** — once a run finishes, the combined, ranked results
  from `top_hits_combined.tsv`, with a relative-strength bar per ligand.
- **Clean / Nuke** — same destructive operations as the CLI flags, each
  gated behind a confirmation dialog before anything is deleted.

### Design notes

- `gui.py` (`PipelineAPI`) never reimplements pipeline logic — it calls the
  exact same `pipeline.py` stage functions the CLI does. See
  [Architecture](#architecture) for the one intentional behavioral
  difference (per-tranche failure isolation).
- The frontend (`gui_assets/index.html`) is a single self-contained file:
  no external network requests, no CDN dependencies, no build step — it's
  loaded straight off disk by pywebview.
- There's no push channel from Python to JS by design; the frontend polls
  `get_status()` on an interval. This keeps the threading model simple (one
  writer thread — the run thread — one reader — the poll loop) at the cost
  of up to ~700ms of display latency, which is a non-issue for a pipeline
  whose stages take seconds to hours.

## Configuration file reference

- **`setup.txt`** — `KEY=VALUE` pairs, one per line (`#` comments and blank
  lines ignored):
  - `recList` — path to `recList.txt`.
  - `geoList` — path to `geoList.txt`.
  - `ligList` — path to `ligList.txt`. Required to exist, but its contents
    aren't consumed directly — per-tranche ligand lists are generated
    automatically in Stage 2 from each tranche's downloaded ligands.
  - `filter_percent` — number in `(0, 100]`; the top slice of ranked
    ligands to keep per tranche.
  - `cores` — integer `>= 1`; degree of parallelism passed as
    `mpirun -np <cores>`.
  - `energy_range` *(optional)* — positive number; passed through to the
    docking binary's config.
- **`tranches.txt`** — TSV with header `curl_script log_p molecular_weight`,
  one row per ZINC tranche. `curl_script` is resolved relative to this
  file's own directory.
- **`recList.txt`** — one receptor `.pdbqt` path per line, resolved
  relative to this file's own directory.
- **`geoList.txt`** — one grid box per line, paired line-for-line with
  `recList.txt`: `center_x center_y center_z size_x size_y size_z`.

## Directory layout produced by a run

```
<workdir>/
├── tranches/
│   └── LP5_MW200/                    # one directory per tranche
│       ├── status.txt                # INIT | DOWNLOADED | UNPACKED | DOCKED | DONE | FAILED_*
│       ├── meta.txt                  # creation time, log_p, molecular_weight
│       ├── inputs/                   # snapshotted curl script + setup.txt
│       ├── download/                 # raw *.pdbqt.gz, as downloaded
│       ├── ligands/                  # decompressed *.pdbqt ligands
│       ├── ligand_list.txt           # paths into ligands/, fed to the docking binary
│       ├── docking/
│       │   └── target1_protein/      # one dir per (receptor, grid box) target
│       │       ├── vinalc.conf
│       │       └── *_out.pdbqt       # docked poses + REMARK VINA RESULT scores
│       ├── logs/                     # stdout/stderr + timing for each subprocess call
│       └── results/
│           ├── ranked_all.tsv        # every docked ligand for this tranche, best-first
│           └── top_hits.tsv          # top filter_percent% for this tranche
└── results/
    └── top_hits_combined.tsv         # top hits across all tranches, merged + re-ranked
```

## Placeholders used for offline testing

Two files in this repo are intentionally placeholders so the pipeline runs
end-to-end without network access or a real docking binary. **Replace both
before running a real screen:**

- **`ZINC-downloader-3D-pdbqt.gz.curl`** — synthesizes two tiny fake
  ligands locally instead of hitting ZINC's servers. Replace with the real
  curl-command file downloaded from the
  [ZINC tranche picker](https://zinc20.docking.org).
- **`protein.pdbqt`** — a minimal placeholder receptor (a few bare CA
  atoms, not a real prepared structure). Replace with a real receptor
  prepared for docking (e.g. via AutoDockTools/MGLTools or Meeko).

You'll also want to point `--vinalc-bin` (CLI) / "Docking binary" (GUI) at
your actual VinaLC (or Vina/smina/QuickVina) install, rather than the
`tests/fixtures/mock_vinalc.py` stand-in used by the automated tests.

## Testing

```bash
python3 -m unittest discover -s tests -v
```

Tests run the full pipeline offline against `tests/fixtures/mock_vinalc.py`,
a stand-in docking binary that writes deterministic fake affinities in the
same output format (`REMARK VINA RESULT:` lines) that Vina/VinaLC produce,
so the real parsing and ranking logic is exercised without needing MPI or a
real docking install. Coverage includes: config parsing/validation edge
cases, a full 5-stage run end-to-end, resumability (re-running a completed
stage is a no-op), and the `--clean`/`--nuke` safety behavior.

## Troubleshooting

- **`FileNotFoundError: recList/geoList/ligList does not exist`** — these
  paths in `setup.txt` are resolved as given (after `~` expansion), not
  relative to `setup.txt`'s own directory. Use absolute paths, or paths
  relative to wherever you run `pipeline.py`/`gui.py` from.
- **`geoList (...) and recList (...) must have the same number of lines`**
  — every receptor in `recList.txt` needs a matching grid box on the same
  line number in `geoList.txt`.
- **Docking stage fails immediately with a "binary not found"-style
  error** — `--vinalc-bin`/`--mpirun-bin` (CLI) or the "Docking binary"/"MPI
  launcher" fields (GUI) must be on `PATH` or given as absolute paths. Use
  `--no-mpirun` (or check "Skip MPI launcher" in the GUI) if you're running
  the docking binary directly without MPI.
- **A tranche is stuck skipping a stage you wanted to re-run** — statuses
  are sticky by design (that's what makes runs resumable). Delete that
  tranche's `status.txt`, or its whole directory under
  `<workdir>/tranches/`, to force it to re-run from scratch; `--clean`
  removes all tranches under a workdir at once.
