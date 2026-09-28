# dockingflow

A resumable virtual-screening pipeline: downloads ZINC ligand tranches,
docks them against a receptor with a parallel AutoDock Vina backend
(VinaLC-style), and ranks/filters the results down to top hits. Comes with
both a CLI (`pipeline.py`) and a desktop GUI (`gui.py`) built on the same
underlying code.

- [Glossary](#glossary)
- [Quick start](#quick-start)
- [Running on the docking server](#running-on-the-docking-server)
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
python3 pipeline.py --map tranches.txt --setup setup.txt --workdir ./run1 \
    --vinalc-bin tests/fixtures/mock_vinalc.py --no-mpirun
```

This uses the placeholder receptor/ligand data and the mock VinaLC checked
into this repo, so it runs fully offline out of the box — see
[Placeholders used for offline testing](#placeholders-used-for-offline-testing)
before pointing it at a real screen.

## Running on the docking server

### With the GUI (recommended; e.g. from MobaXterm)

1. SSH into the server with MobaXterm (X11-Forwarding is on by default in
   MobaXterm SSH sessions), `cd` into this repo, and run:
   ```bash
   bash gui.sh
   ```
   A DockingFlow window opens on your own screen. Everything in it —
   machine detection, file browsing, downloads, docking — happens on the
   server.
2. In the window, top to bottom:
   - **Import ZINC downloader…** — pick the `.curl` file from
     CartBlanche22 (copy it to the server first, e.g. by dragging it into
     MobaXterm's file browser). It's split into one download script per
     tranche and the Tranches map field is filled in.
   - **Docking targets** — Browse to your prepared receptor `.pdbqt`, and
     enter the grid box center and size (Å) around the binding pocket.
     Add more receptors with **+ Add receptor**.
   - **Docking binary** — `vinalc` if it's on PATH, otherwise its full
     path.
   - **Resources** — the server is detected automatically and the sliders
     are pre-set to a recommendation; adjust if you're sharing the machine.
   - **Validate**, then **Run pipeline**. Progress, the log, and the top
     hits update live on the right.
3. Close the window or MobaXterm whenever you like: the run keeps going on
   the server. Run `bash gui.sh` again later to reopen the window and see
   live progress. `bash gui.sh stop` stops the GUI server (and any run in
   progress — re-running resumes where it stopped).

The GUI keeps its settings in `project/` (created on first launch from the
repo's example files, and ignored by git), so `git pull` for code updates
never conflicts with your settings.

How it works: `gui.sh` starts `gui.py --web` detached in the background (it
owns the run; log in `gui_server.log`) and then opens a window onto it
through MobaXterm's X server. That window is a Tkinter app
([`gui_tk.py`](gui_tk.py)), which is fast over X11 because Tk sends small
drawing commands rather than rendered pixels; it needs Python's `tkinter`
(the `python3-tk` package on some Linux installs — `server_check.sh`
reports whether it's there). Without tkinter, `gui.sh` falls back to a
browser-based window — much slower over X11 — using a browser on the
server, or pywebview installed by `bash gui.sh setup`.
`bash gui.sh browser` forces the browser-based window. If no window can be opened, `gui.sh` prints how to
reach the GUI from your own browser through a MobaXterm SSH tunnel
instead. The web server listens on localhost only and requires the random
token in its URL, so other users on a shared server can't control it.

### With the command line

1. **Check the machine.** Copy the repo over and run
   `bash server_check.sh /path/to/workdir 2>&1 | tee server_report.txt`.
   It reports CPU cores (physical vs. hyperthreads), RAM, current load,
   free disk space and inodes, Python, `vinalc`/`mpirun` (and an
   `mpirun -np 2` smoke test), job schedulers, and whether ZINC is
   reachable. It's read-only. Its last block (`DF_*=` lines) is what the
   GUI's **Load server report** button reads.
2. **Split the ZINC22 downloader into tranches.** CartBlanche22 gives one
   file covering every tranche you picked:
   ```bash
   python3 zinc_split.py ZINC22-downloader-3D-pdbqt.tgz.curl --out-dir zinc22_scripts --map tranches.txt
   ```
   This writes one download script per tranche (e.g.
   `zinc22_scripts/H04M000.curl`) and a `curl_script tranche` mapping
   file.
3. **Pick a CPU/memory budget** — either in the GUI's Resources panel (then
   **Save to setup file**), or by editing `cores=` / `memory_gb=` in
   `setup.txt` directly. See [CPU and memory](#cpu-and-memory).
4. **Replace the placeholders** — a real prepared receptor in
   `recList.txt`, and its grid box in `geoList.txt`.
5. **Run**, detached so it survives logging out:
   ```bash
   nohup python3 pipeline.py --map tranches.txt --setup setup.txt --workdir /data/run1 \
       --vinalc-bin /path/to/vinalc > run1.log 2>&1 &
   ```
   Re-running the same command resumes where it stopped.

### CPU and memory

VinaLC runs one master rank (hands out ligands) plus worker ranks, and each
worker runs about `exhaustiveness` (default 8) search threads of its own.
So `mpirun -np <cores>` would oversubscribe the machine ~8-fold. Instead
the pipeline launches `cores // exhaustiveness` workers + 1 master (e.g.
`cores=192` → `-np 25`).

Each worker also builds a double-precision grid map per ligand atom type
spanning the whole search box, so memory per worker grows with the box
volume: roughly 0.8 GB for an 80 Å cube at the default 0.375 Å granularity
(boxes of ~20–30 Å, typical for a known pocket, need far less). If
`memory_gb` is set, the worker count is capped to fit that budget too.
`mpi_ranks=` in setup.txt overrides both.

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
   snapshotted curl script, then verifies at least one ligand file
   (`*.pdbqt.tgz`, `*.pdbqt.gz` or `*.pdbqt`) was produced. Logs
   stdout/stderr to `<tranche>/logs/download.*.log`. If some commands fail
   but others succeed (ZINC regularly 404s on individual archives), the
   tranche continues and the failures are listed in
   `<tranche>/download_failures.txt`.
3. **Stage 2 — unpack** (`unpack_tranche`): reads every molecule out of
   the downloads (ZINC22 `.pdbqt.tgz` archives hold one small pdbqt per
   molecule) and writes them into `<tranche>/ligands/ligands_00001.pdbqt`
   etc., 10,000 per file, **each wrapped in `MODEL`/`ENDMDL`**: VinaLC
   treats each such block as one docking job and ignores anything outside
   them. Also writes `ligand_list.txt` (VinaLC's `--ligList`) and
   `ligand_index.tsv`, mapping VinaLC's `LIGAND <n>` numbering back to
   ZINC ids. Unreadable downloads are skipped and listed in
   `unpack_warnings.txt`.
4. **Stage 3 — dock** (`dock_tranche`): one VinaLC run per tranche, from
   `<tranche>/docking/`:
   `mpirun -np <ranks> vinalc --recList recList.txt --ligList ligList.txt --geoList geoList.txt --exhaustiveness ... --num_modes ... --energy_range ...`.
   VinaLC itself loops over every receptor × ligand. The exact command is
   saved to `docking/command.txt`; results land in
   `docking/recList.txt_ligList.txt.pdbqt.gz` (+ `.log.gz`).
5. **Stage 4 — rank + filter** (`rank_and_filter_tranche`): streams
   VinaLC's combined output (records arrive in completion order, labelled
   only `LIGAND <n>`), takes each record's first `REMARK VINA RESULT:`
   (poses are best-first), maps it back to a ZINC id, keeps each ligand's
   best score across receptors, and writes `results/ranked_all.tsv`,
   `results/top_hits.tsv` (top `filter_percent`%, at least 1 ligand) and
   `results/summary.txt` (submitted vs. successfully docked).
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
| `--vinalc-bin PATH` | VinaLC binary (default `vinalc`, looked up on PATH; relative paths are fine). |
| `--mpirun-bin PATH` | MPI launcher to wrap VinaLC with (default `mpirun`). |
| `--no-mpirun` | Invoke the docking binary directly. Only for the test mock: real VinaLC exits unless it has at least 2 MPI ranks. |
| `--clean` | Delete generated tranche outputs under the workdir; keeps the workdir itself. Does not require `--map`/`--setup` to be valid. |
| `--nuke` | Delete the entire workdir. **Destructive and irreversible.** Does not require `--map`/`--setup` to be valid. |

## GUI

A desktop app (not a browser tab) built with [pywebview](https://pywebview.flowrl.com/):
Python hosts the pipeline logic and exposes it to a bundled HTML/CSS/JS page
rendered in a native window.

### Running

- **On the docking server:** `bash gui.sh` — see
  [Running on the docking server](#with-the-gui-recommended-eg-from-mobaxterm).
  Needs nothing beyond Python's standard library on the server (plus a
  browser, or `bash gui.sh setup`, for the window).
- **As a local desktop app** (runs the pipeline on this machine):
  ```bash
  python3 -m venv .venv
  .venv/bin/pip install pywebview
  .venv/bin/python3 gui.py
  ```
- **Web mode by hand:** `python3 gui.py --web [--port 8765]` prints an SSH
  tunnel command and a `http://localhost:8765/?token=...` URL.

### What it does

- **Configuration panel** — set the setup file, tranches map, work
  directory, docking binary, and MPI launcher, with native file/folder
  pickers. Opens pre-filled with this repo's own `setup.txt`/`tranches.txt`.
- **Resources panel** — shows the machine's physical/logical cores, RAM
  and load, either **detected** (the machine the GUI runs on, done
  automatically at startup) or from a **server report** (the output of
  `server_check.sh`, for sizing a remote server from your laptop). Sliders
  pick how many cores and how much RAM the run may use, pre-set to a
  recommendation (physical cores minus current load and a ~5% reserve; 80%
  of available RAM), and show live what that means: the `mpirun -np`,
  estimated memory, and whether cores or memory is the limit. The budget
  applies to Validate/Run from the GUI; **Save to setup file** writes it
  as `cores=`/`memory_gb=` for CLI runs on the server.
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
  no external network requests, no CDN dependencies, no build step. The
  same page runs in both modes: in the desktop app it calls
  `pywebview.api`; served by `web_server.py` it gets `window.DF_WEB = true`
  and routes the same calls over HTTP (`POST /api/<method>`), with an
  in-page file picker (`PipelineAPI.list_dir`) replacing native dialogs.
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
  - `ligList` *(optional, unused)* — per-tranche ligand lists are
    generated automatically in Stage 2.
  - `filter_percent` — number in `(0, 100]`; the top slice of ranked
    ligands to keep per tranche.
  - `cores` — integer `>= 1`; CPU cores the run may use. See
    [CPU and memory](#cpu-and-memory) for how this becomes `mpirun -np`.
  - `memory_gb` *(optional)* — RAM budget; caps the number of VinaLC
    workers.
  - `energy_range` *(optional, default 3)*, `exhaustiveness` *(8)*,
    `num_modes` *(9)*, `granularity` *(0.375)*, `seed` — passed to VinaLC
    as the matching `--flags`.
  - `mpi_ranks` *(optional)* — explicit `mpirun -np`, overriding
    `cores`/`memory_gb`.

  Relative paths are resolved against setup.txt's own directory.
- **`tranches.txt`** — whitespace-separated, one row per tranche, with
  header either `curl_script tranche` (ZINC22 codes like `H04M000`, as
  written by `zinc_split.py`) or `curl_script log_p molecular_weight`
  (older ZINC20-style). `curl_script` is resolved relative to this file's
  own directory.
- **`recList.txt`** — one receptor `.pdbqt` path per line, resolved
  relative to this file's own directory.
- **`geoList.txt`** — one grid box per line, paired line-for-line with
  `recList.txt`: `center_x center_y center_z size_x size_y size_z`.

## Directory layout produced by a run

```
<workdir>/
├── tranches/
│   └── H04M000/                      # one directory per tranche
│       ├── status.txt                # INIT | DOWNLOADED | UNPACKED | DOCKED | DONE | FAILED_*
│       ├── meta.txt                  # creation time, tranche, log_p, size bin
│       ├── inputs/                   # snapshotted curl script + setup.txt
│       ├── download/                 # raw *.pdbqt.tgz, as downloaded
│       ├── download_failures.txt     # download commands that failed (if any)
│       ├── ligands/                  # ligands_00001.pdbqt ...: MODEL-wrapped molecules
│       ├── ligand_list.txt           # paths into ligands/ (VinaLC --ligList)
│       ├── ligand_index.tsv          # VinaLC "LIGAND <n>" -> ZINC id
│       ├── unpack_warnings.txt       # unreadable downloads (if any)
│       ├── docking/
│       │   ├── recList.txt, geoList.txt, ligList.txt, command.txt
│       │   ├── recList.txt_ligList.txt.pdbqt.gz   # all poses + REMARK VINA RESULT scores
│       │   └── recList.txt_ligList.txt.log.gz
│       ├── logs/                     # stdout/stderr + timing for each subprocess call
│       └── results/
│           ├── ranked_all.tsv        # ligand, affinity, receptor — best-first
│           ├── top_hits.tsv          # top filter_percent% for this tranche
│           └── summary.txt           # submitted vs. docked counts
└── results/
    └── top_hits_combined.tsv         # ligand, affinity, receptor, tranche — merged + re-ranked
```

## Placeholders used for offline testing

Two files in this repo are intentionally placeholders so the pipeline runs
end-to-end without network access or a real docking binary. **Replace both
before running a real screen:**

- **`ZINC-downloader-3D-pdbqt.gz.curl`** — synthesizes a tiny archive in
  the real ZINC22 `.pdbqt.tgz` layout locally instead of hitting ZINC's
  servers. For a real screen, run `zinc_split.py` on the downloader from
  [CartBlanche22](https://cartblanche22.docking.org) and use the
  `tranches.txt` it writes.
- **`protein.pdbqt`** — a minimal placeholder receptor (a few bare CA
  atoms, not a real prepared structure). Replace with a real receptor
  prepared for docking (e.g. via AutoDockTools/MGLTools or Meeko).

You'll also want to point `--vinalc-bin` (CLI) / "Docking binary" (GUI) at
your actual VinaLC install, rather than the `tests/fixtures/mock_vinalc.py`
stand-in used by the automated tests. The pipeline speaks VinaLC's command
line and output format specifically; plain Vina/smina/QuickVina won't work
without changing `build_vinalc_command` and `parse_vinalc_poses`.

## Testing

```bash
python3 -m unittest discover -s tests -v
```

Tests run the full pipeline offline against `tests/fixtures/mock_vinalc.py`,
a stand-in that takes VinaLC's real flags and writes VinaLC's real output
file and record format (out of order, `LIGAND <n>` labels) with
deterministic fake affinities, so the real parsing and ranking logic is
exercised without MPI or a real docking install.
`tests/fixtures/H04M000-N-aaaaaa.pdbqt.tgz` is a real 4-molecule ZINC22
archive. Coverage includes config parsing/validation, a full run
end-to-end, resumability, partial download failures, ZINC22/ZINC20 unpack
formats, `zinc_split.py`, CPU/memory planning and server-report parsing,
and the `--clean`/`--nuke` safety behavior.

## Troubleshooting

- **`FileNotFoundError: recList/geoList does not exist`** — relative
  paths in `setup.txt` are resolved against `setup.txt`'s own directory
  (not the directory you launch from).
- **Docking fails with `Error: Total process less than 2`** — VinaLC needs
  at least 2 MPI ranks; don't use `--no-mpirun` with the real binary.
- **`mpirun` refuses to start as root / "not enough slots"** — Open MPI
  needs `--allow-run-as-root` as root, and may cap `-np` at the core count
  it detects. `server_check.sh` runs an `mpirun -np 2` smoke test to catch
  launcher problems before a real run.
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
