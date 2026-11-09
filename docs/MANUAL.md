# DockingFlow — Lab Manual

A practical guide to running structure-based virtual screens with DockingFlow
on the lab's docking server. No prior experience with the tool is assumed.

If you only need to get something running today, read
[2. Your first screen](#2-your-first-screen) and stop there.

**Contents**

1. [What DockingFlow does](#1-what-dockingflow-does)
2. [Your first screen (about 30 minutes)](#2-your-first-screen-about-30-minutes)
3. [Preparing a receptor](#3-preparing-a-receptor)
4. [The window, tab by tab](#4-the-window-tab-by-tab)
5. [Running, stopping and resuming](#5-running-stopping-and-resuming)
6. [Results](#6-results)
7. [Planning a real campaign](#7-planning-a-real-campaign)
8. [Command line](#8-command-line)
9. [Settings reference (`setup.txt`)](#9-settings-reference-setuptxt)
10. [Troubleshooting](#10-troubleshooting)
11. [Glossary](#11-glossary)
12. [How it works inside](#12-how-it-works-inside)
13. [Known limitations](#13-known-limitations)

---

## 1. What DockingFlow does

Virtual screening asks a computer which molecules out of a large catalogue are
most likely to bind a target protein. DockingFlow automates the whole loop:

```
   choose compounds  ->  download  ->  prepare  ->  dock  ->  rank  ->  inspect hits
     (ZINC22 or        (per tranche,   (3D, atom   (VinaLC    (best     (poses for
      ChEMBL)           resumable)      types)      on N cores) score)    PyMOL)
```

The docking engine is **VinaLC**, a version of AutoDock Vina that spreads work
across many CPU cores at once. DockingFlow does everything around it: picking
compounds, fetching them, converting formats, sizing the job to the machine,
restarting cleanly after interruptions, and collecting the results.

**Key ideas**

- A **screen** is split into **tranches** — batches of a few thousand to a few
  hundred thousand compounds. Each tranche is downloaded, prepared, docked and
  ranked in turn, and each records its progress on disk. If a run is
  interrupted, finished tranches are never redone.
- Everything happens **on the server**. The window on your screen is only a
  remote control; you can close it and the run carries on.
- A **work directory** holds everything a screen produces. It can get large
  (roughly 1–5 GB per million compounds), so put it on a disk with room.

---

## 2. Your first screen (about 30 minutes)

This runs a deliberately tiny screen end to end, so you learn the tool before
committing to a long job.

### 2.1 Connect and open the window

1. Open **MobaXterm** and start your usual SSH session to the docking server.
   (X11 forwarding is on by default in MobaXterm SSH sessions; it is what lets
   the window appear on your PC.)
2. Go to the DockingFlow folder and launch it:

   ```bash
   cd ~/dockingflow        # wherever the repository lives
   bash gui.sh
   ```

A window titled **DockingFlow** opens on your screen. If it does not, see
[Troubleshooting](#10-troubleshooting) — `gui.sh` prints what to fix.

![The DockingFlow window just after opening](images/01-overview.png)

The window has two halves: the numbered tabs on the left are what you set up,
and the panels on the right show what is happening.

### 2.2 Work through the tabs

Click **How to use** (top right) for a short tour, and the **?** next to any
setting for an explanation of that setting.

**Tab 1 — Ligands.** Leave the source as **ZINC22** and the library as
**Quick test**. That is one small tranche of about a thousand compounds. Check
the **Work directory** points somewhere with free space.

**Tab 2 — Targets.** Click **Browse** and pick your prepared receptor
(`.pdbqt` — see [section 3](#3-preparing-a-receptor) if you do not have one).
Then set the search box:

- Know where the ligand binds? **Center on known site…** and pick a structure
  file containing a ligand sitting in that pocket.
- Don't know? **Fit to whole protein** (blind docking).

**Tab 3 — Docking.** Check the docking program is `vinalc` (or give its full
path). Leave the other settings alone for now.

**Tab 4 — Resources.** The server is measured automatically and the sliders
are pre-set to a sensible share of it. Lower them if colleagues are using the
machine.

### 2.3 Validate, then run

Click **Validate**. It checks every input and reports what the run will do.
Fix anything it complains about, then click **Run pipeline**.

The right-hand side now fills in: tranche progress at the top, a live log in
the middle, and hits at the bottom as they are found.

![A run in progress](images/02-running.png)

When it finishes, the strongest predicted binders are listed in **Top hits**,
and the full results are in your work directory under `results/`.

### 2.4 Look at a hit

Click **Export poses…** on the Top hits panel and accept the suggested number.
DockingFlow writes each hit's best docked pose to `results/poses/`, along with
a copy of the receptor. Copy that folder to your PC (drag it out of MobaXterm's
file browser on the left) and open it in PyMOL or ChimeraX to see how the
compounds sit in the pocket.

That is the whole workflow. Everything else in this manual is detail.

---

## 3. Preparing a receptor

DockingFlow does **not** prepare receptors; you do that once, beforehand.

1. **Get a structure.** Download your target from the
   [Protein Data Bank](https://www.rcsb.org), or use your own model.
2. **Clean it.** Remove waters, ions and any co-crystallised ligand (keep the
   ligand in a separate file — it is the easiest way to define the search box).
   Add missing side chains if needed.
3. **Convert to PDBQT**, which adds hydrogens, charges and AutoDock atom types.
   Any of these works:

   ```bash
   # ADFR suite (recommended)
   prepare_receptor -r protein.pdb -o protein.pdbqt

   # MGLTools / AutoDockTools
   prepare_receptor4.py -r protein.pdb -o protein.pdbqt -A checkhydrogens

   # Meeko
   mk_prepare_receptor.py -i protein.pdb -o protein.pdbqt
   ```

**A PDBQT receptor is not the same as a PDBQT ligand.** VinaLC's receptor
reader accepts only `ATOM`/`HETATM`, `REMARK` and `TER` lines. A file with
`ROOT`, `BRANCH` or `TORSDOF` lines is in *ligand* format and every docking job
will fail. DockingFlow checks this when you click Validate and names the
offending line, but it is easier to get right the first time.

### Choosing the search box

| Situation | What to do | Typical size |
|---|---|---|
| Binding site known (co-crystal ligand, or a published site) | **Center on known site…**, pick the ligand file | 20–30 Å per side |
| Site unknown | **Fit to whole protein** | up to 80 Å per side |

Smaller boxes dock faster and more accurately. For a blind-docking box, raise
**exhaustiveness** to 16–32 in tab 3 so the larger space is searched properly.

---

## 4. The window, tab by tab

### 4.1 Tab 1 — Ligands

Where your compounds come from and how many of them.

```
┌─ 1 Ligands ─────────────────────────────┐
│ Compound source            [?]          │
│ [ ZINC22 (ready-to-dock 3D)        v ]  │
│                                         │
│ ZINC22 library             [?]          │
│ [ Quick test (~1,000 molecules)    v ]  │
│ Heavy atoms  [10 v] to [10 v]      [?]  │
│ logP         [3.1 v] to [3.1 v]    [?]  │
│ Charges      [ Neutral only        v ]  │
│ [x] 0  [ ] -1  [ ] +1  ...              │
│                                         │
│ ┌─────────────────────────────────────┐ │
│ │ ~1,001 molecules                    │ │
│ │ in 1 tranche(s).                    │ │
│ │ [ Create tranche list ]             │ │
│ └─────────────────────────────────────┘ │
│                                         │
│ Tranche list               [?]          │
│ [ .../project/zinc22_tranches.txt ][…]  │
│ Work directory             [?]          │
└─────────────────────────────────────────┘
```

![The Ligands tab with ZINC22 selected](images/03-ligands-zinc.png)

**Compound source** — the one choice that changes everything else:

| | **ZINC22** | **ChEMBL** |
|---|---|---|
| What it is | Purchasable compounds (billions) | Bioactive compounds with measured activity data (~1.9 million small molecules) |
| Format supplied | 3D, ready to dock | 2D only |
| Extra software | none | **Open Babel** on the server |
| Speed | fastest route to docking | slower: every compound must be built in 3D first |
| Use it when | exploring purchasable chemical space | you want known-bioactive compounds, or to match a published ChEMBL screen |

**ZINC22 presets**

| Preset | Heavy atoms | logP | Roughly |
|---|---|---|---|
| Quick test | 10 | 3.1 | 1 thousand |
| Fragments (Rule of Three) | 8–19 | −1 to 3 | hundreds of millions |
| Lead-like | 17–25 | −1 to 3.5 | billions |
| Drug-like (Lipinski) | 17–29 | −1 to 5 | billions |

ZINC22 sizes molecules by **heavy-atom count** (non-hydrogen atoms) rather than
molecular weight; multiply by about 14 for an approximate weight.

**ChEMBL presets**

| Preset | Weight (Da) | logP | Roughly |
|---|---|---|---|
| Quick test | 200–201 | — | 900 |
| Fragments (Rule of Three) | ≤ 300 | ≤ 3 | 230 thousand |
| Lead-like | 250–350 | ≤ 3.5 | 370 thousand |
| Drug-like (Lipinski) | ≤ 500 | ≤ 5 | over a million |
| Whole library | — | — | 1.9 million |

![The Ligands tab with ChEMBL selected](images/04-ligands-chembl.png)

With ChEMBL you also choose **how** it is fetched:

- **Filtered selection (resumable)** — compounds come from ChEMBL's web API in
  pages of 1000, about 10 seconds per page, split into one resumable tranche
  per 25 Da weight bin. Use this unless you want the entire database.
- **Whole release in one file** — downloads the complete release (about 1 GB)
  from EBI. Far quicker for all 1.9 million compounds, but it arrives as a
  *single* tranche: no intermediate checkpoints, your filters are ignored, and
  docking runs as one very long job.

**Check Open Babel** (ChEMBL only) converts two sample compounds and reports
what your server's Open Babel actually did — its version, whether 3D
coordinates came out, whether compound names survived. **Run this once on a new
server before a ChEMBL screen.**

The count under the filters is live, and **Create tranche list** writes the
download plan. Nothing is downloaded until you run the pipeline.

> **Importing from the ZINC website.** The *Import a CartBlanche22 file*
> button accepts a download script you got from
> [cartblanche22.docking.org](https://cartblanche22.docking.org). DockingFlow
> repairs two problems in those files automatically: ZINC's own script gives
> every archive in a tranche the same output name (so downloads overwrite each
> other), and its links often fail. See
> [the README](../README.md#zinc22-downloads) for details. Using the built-in
> picker instead avoids the issue entirely.

### 4.2 Tab 2 — Targets

One or more receptors, each with its own search box.

```
┌─ 2 Targets ─────────────────────────────┐
│ Receptors and grid boxes   [?]          │
│ ┌─ Receptor 1 ────────────────────────┐ │
│ │ [ /path/protein.pdbqt ][Browse][Rm] │ │
│ │            x        y        z      │ │
│ │ Center  [ 97.4 ] [127.5 ] [134.0 ]  │ │
│ │ Size(Å) [ 22.0 ] [ 22.0 ] [ 22.0 ]  │ │
│ │ [Fit to whole protein][?]           │ │
│ │ [Center on known site...][?]        │ │
│ └─────────────────────────────────────┘ │
│ [ + Add receptor ]  [ Save targets ]    │
└─────────────────────────────────────────┘
```

![The Targets tab](images/05-targets.png)

Add more receptors to dock the same compounds against several structures (for
example different conformations, or related proteins). Each hit records which
receptor gave its best score.

Box sides are capped at 80 Å, which is enough to cover a typical protein.

### 4.3 Tab 3 — Docking

The docking program and its settings. The defaults are sensible; change them
only with a reason.

| Setting | Default | What it does |
|---|---|---|
| Keep top % | 10 | How much of the ranking is kept as "hits". Everything is still saved. |
| Keep top % of | whole screen | Whether that percentage is of the whole screen or of each tranche separately. |
| Exhaustiveness | 8 | Search effort per compound. Higher is more reliable and proportionally slower. Use 16–32 for blind docking. |
| Poses per ligand | 9 | How many poses are reported. Only the best is ranked. |
| Energy range | 3 | kcal/mol window of reported poses. |
| Random seed | blank | Set a number to make a run reproducible. |
| Grid spacing | 0.375 | Vina's grid resolution. Smaller is slower and uses much more memory. |

For ChEMBL screens this tab also holds the ligand-preparation settings
(Open Babel path, 3D generation effort, minimisation steps).

### 4.4 Tab 4 — Resources

How much of the server this run may use.

![The Resources tab](images/06-resources.png)

The panel shows the detected machine and what your choice means in practice:

```
VinaLC will run mpirun -np 25: 1 coordinator + 24 workers x ~8 threads.
Estimated memory: 20.3 GB (~0.8 GB per worker for your grid box).
Limited by: cores.
```

Why not one worker per core? Each VinaLC worker itself runs about
*exhaustiveness* threads, so one worker per core would oversubscribe the
machine roughly eightfold. DockingFlow works out the right number for you.

**Load server report…** lets you size a run for a *different* machine, using
the output of `bash server_check.sh` from that machine.

---

## 5. Running, stopping and resuming

**Validate** first: it checks every input without downloading or docking
anything, and saves your settings. Then **Run pipeline**.

Each tranche moves through these states:

```
Waiting -> Downloaded -> Unpacked -> Docked -> Done
```

(ChEMBL screens spend most of "Unpacked" building 3D structures.)

- **Closing the window does not stop the run.** Neither does closing
  MobaXterm, or losing your connection. Reopen with `bash gui.sh` to check on
  it; the tables repopulate from the work directory.
- **Stop** ends the current step and all of its processes. Everything already
  finished is kept, and the button becomes **Resume**. The interrupted step is
  redone from its start — so stopping mid-docking means that tranche re-docks,
  while finished tranches are skipped.
- **A failed tranche does not stop the screen.** It is recorded and the run
  moves on. Select it and click **Details…** to see why.
- `bash gui.sh stop` shuts down the server *and* any run. Use it when you are
  finished for the day and want the machine free.

### Watching a long run

Select any tranche and click **Details…** (or double-click it) for a single
report: how many compounds were downloaded, converted, docked and ranked, the
exact VinaLC command, any errors VinaLC printed, and the tails of its logs.
**Copy report** puts it on the clipboard — the quickest way to ask someone for
help.

![The per-tranche details report](images/07-details.png)

---

## 6. Results

### Where everything lands

```
<work directory>/
├── results/
│   ├── top_hits_combined.tsv     <- the hits: your main result
│   ├── ranked_all_combined.tsv   <- every docked compound, best first
│   └── poses/                    <- after "Export poses"
│       ├── 001_CHEMBL123.pdbqt   <- each hit's best pose
│       ├── top_hits_poses.pdbqt  <- all hits in one file
│       └── receptor_protein.pdbqt
└── tranches/<tranche>/
    ├── status.txt                <- how far this tranche got
    ├── download/                 <- raw downloads
    ├── ligands/                  <- prepared compounds fed to VinaLC
    ├── ligand_index.tsv          <- VinaLC's "LIGAND n" -> real compound id
    ├── docking/                  <- VinaLC's config, command and raw output
    ├── logs/                     <- stdout/stderr of every command
    └── results/
        ├── ranked_all.tsv        <- this tranche, best first
        └── top_hits.tsv
```

`top_hits_combined.tsv` is a tab-separated table you can open in Excel or read
with pandas:

| ligand | affinity_kcal_mol | receptor | tranche |
|---|---|---|---|
| CHEMBL1234 | −11.2 | protein | MW300-325 |

### Reading the scores

Affinity is in **kcal/mol, and more negative means stronger predicted
binding**. Treat these as a *ranking*, not a measurement: Vina-family scores
are good at sorting compounds but poor as absolute binding energies, and the
error on any single number is roughly 2–3 kcal/mol.

Practical advice:

- **Do not simply take the top 10.** Look at tens to hundreds of hits, check
  the poses make chemical sense, and watch for recurring scaffolds.
- **Large molecules score better** almost automatically, because more atoms
  make more contacts. Compare within a narrow size range, or use ligand
  efficiency (affinity ÷ heavy atoms).
- **Check the pose, not just the number.** Export poses and look at whether
  the compound sits in the real pocket making sensible interactions.
- **Expect false positives.** Virtual screening is a filter that turns
  millions of compounds into a testable shortlist, not a prediction of
  activity.

### Hits as they appear

The Top hits table updates after each tranche finishes, so you can start
looking at early results while a long screen continues.

---

## 7. Planning a real campaign

### How long will it take?

Docking dominates. A rough estimate:

```
core-seconds  ≈  compounds × exhaustiveness/8 × (box volume factor) × ~10
wall-clock    ≈  core-seconds ÷ (cores you allocated)
```

On 192 cores, with default settings and a 22 Å box, expect very roughly:

| Compounds | Rough wall-clock |
|---|---|
| 1 thousand | minutes |
| 100 thousand | hours |
| 1 million | days |
| 10 million+ | weeks |

Add to that, for ChEMBL only, 3D conversion at roughly 0.1–1 s per compound
per core, and for the ChEMBL API route about 10 s per 1000 compounds of
downloading.

**Always measure first.** Run one small tranche, look at the timing in
`docking/dock_timing.txt`, and scale up from there. That single number beats
any estimate in this manual.

### Disk

Budget roughly 1–5 GB per million compounds, more if you keep raw downloads.
Check free space with `bash server_check.sh /path/to/workdir` before a big run.

### A sensible campaign

1. **Pilot** — Quick test preset, confirm the whole chain works end to end.
2. **Calibrate** — one real tranche; measure speed, check hits look sane.
3. **Scale** — widen the filters to the library you actually want.
4. **Triage** — inspect poses of the top hits, cluster by scaffold, check
   availability and drug-likeness.
5. **Rescore (optional)** — re-dock the shortlist with higher exhaustiveness,
   or with a second program, before ordering anything.

### Being a good neighbour

The server is shared. Leave cores free (the Resources tab already reserves
some), and tell colleagues before starting a multi-day run.

---

## 8. Command line

Everything the GUI does is available without it — useful for scripting and for
reproducing a run exactly.

```bash
python3 pipeline.py --map project/zinc22_tranches.txt \
                    --setup project/setup.txt \
                    --workdir /data/screens/run1 \
                    --vinalc-bin /path/to/vinalc
```

Run it detached so it survives logging out:

```bash
nohup python3 pipeline.py --map … --setup … --workdir … > run1.log 2>&1 &
```

| Flag | Meaning |
|---|---|
| `--map PATH` | Tranche list (required) |
| `--setup PATH` | Settings file (required) |
| `--workdir PATH` | Output directory (required) |
| `--vinalc-bin PATH` | VinaLC binary (default `vinalc`) |
| `--mpirun-bin PATH` | MPI launcher (default `mpirun`) |
| `--no-mpirun` | Run the binary directly — testing only; real VinaLC needs ≥ 2 MPI ranks |
| `--only-stage0` / `--only-download` / `--only-unpack` | Stop after that stage |
| `--clean` | Delete generated tranche outputs, keep the work directory |
| `--nuke` | Delete the whole work directory (irreversible) |

Ctrl-C (or `kill`) stops cleanly; re-running the same command resumes.

Building a tranche list without the GUI:

```bash
# ZINC22, from a CartBlanche download script
python3 zinc_split.py ZINC22-downloader-3D-pdbqt.tgz.curl \
        --out-dir zinc22_scripts --map tranches.txt

# ChEMBL, from Python
python3 -c "import chembl, pathlib; \
  print(chembl.build_tranche_map(pathlib.Path('chembl_scripts'), \
        pathlib.Path('tranches.txt'), mw_min=250, mw_max=350, logp_max=3.5))"
```

Other useful commands:

```bash
bash gui.sh            # open the window (starts the server if needed)
bash gui.sh status     # is the server running? print its URL
bash gui.sh restart    # after a git pull, so the server picks up new code
bash gui.sh stop       # stop the server and any run
bash gui.sh browser    # use a browser window instead of the Tk one
bash server_check.sh   # report the machine's cores, memory, disk and software
```

---

## 9. Settings reference (`setup.txt`)

The GUI writes this file for you; you rarely need to edit it. It is
`KEY=VALUE`, one per line, `#` for comments. Relative paths resolve against the
file's own directory.

**Required**

| Key | Meaning |
|---|---|
| `recList` | File listing receptor `.pdbqt` paths, one per line |
| `geoList` | File listing search boxes, one per line, paired with `recList` |
| `filter_percent` | Percentage of the ranking kept as hits (0–100] |
| `cores` | CPU cores this run may use |

**Docking**

| Key | Default | Meaning |
|---|---|---|
| `rank_scope` | `overall` | `overall` or `per_tranche` |
| `exhaustiveness` | 8 | Search effort per compound |
| `num_modes` | 9 | Poses reported per compound |
| `energy_range` | 3 | kcal/mol window of reported poses |
| `granularity` | 0.375 | Grid spacing, Å |
| `seed` | — | Fixed random seed |
| `memory_gb` | — | Memory budget; caps the number of workers |
| `mpi_ranks` | — | Explicit `mpirun -np`, overriding `cores`/`memory_gb` |

**Ligand preparation (SDF sources such as ChEMBL only)**

| Key | Default | Meaning |
|---|---|---|
| `obabel_bin` | `obabel` | Open Babel binary |
| `gen3d` | `fast` | 3D effort: `fastest`, `fast`, `med`, `slow`, `slowest` |
| `forcefield` | `MMFF94` | Minimisation forcefield |
| `minimize_steps` | 0 | Minimisation steps; 0 skips it |
| `prep_ph` | — | Protonate for this pH |
| `prep_jobs` | cores | Parallel Open Babel processes |

---

## 10. Troubleshooting

**The window will not open.**
`gui.sh` prints the reason. Usually one of:
- *X11 forwarding is off* — in MobaXterm, edit the session: SSH → Advanced SSH
  settings → tick X11-Forwarding, and make sure the X server button (top
  right) is on. Reconnect.
- *tkinter is missing* — ask your administrator for the `python3-tk` package.
  Until then `gui.sh` falls back to a slower browser window.
- *Nothing works* — `gui.sh` prints instructions for reaching the GUI from
  your own browser through a MobaXterm tunnel.

**"unknown API method: …"**
The server is running code from before the last `git pull`. Run
`bash gui.sh restart`. The window also detects this and offers to restart.

**Every tranche fails: "none of the N docking job(s) produced a pose".**
VinaLC exits normally even when every job fails, so the real reason is in its
error output — the message quotes it, and **Details…** shows more. The usual
cause is the receptor being in ligand format (`ROOT`/`TORSDOF` lines) or a
plain `.pdb` renamed to `.pdbqt`. See [section 3](#3-preparing-a-receptor).

**"Error: Total process less than 2".**
Real VinaLC needs at least 2 MPI ranks. Do not use `--no-mpirun` with it.

**Downloads fail for some compounds.**
Normal. ZINC lists archives that no longer exist; the failures are recorded in
`download_failures.txt` and the run continues. If *every* file fails, check the
server can reach the internet (`bash server_check.sh`).

**"ZINC's tranche server didn't respond".**
CartBlanche has regular outages. Wait and retry, or switch the source to
ChEMBL.

**"Open Babel is not installed".**
Only ChEMBL screens need it. Install it
(`conda install -c conda-forge openbabel`) or ask your administrator for the
`openbabel` package, then click **Check Open Babel**.

**The window freezes or stops responding.**
It should not: every call that reaches ZINC or ChEMBL, or runs Open Babel,
happens on a background thread while the window stays live (you will see
"Counting..." or "the window stays usable"). If it does freeze, note what you
clicked and report it — and check you are on the current code with
`git pull` followed by `bash gui.sh restart`, since an older window blocked
on those calls.

**The run is slower than expected.**
Check the Resources tab's plan line. Common causes: too few cores allocated, a
very large (blind-docking) box, or high exhaustiveness.

**I need to start completely over.**
**Clean outputs** deletes the tranche results but keeps the directory;
**Nuke** deletes the whole work directory. Both ask first, and both are
irreversible.

---

## 11. Glossary

**Affinity** — the predicted strength of binding, in kcal/mol; more negative is
stronger.

**Docking** — predicting how a small molecule sits in a protein's binding site,
and scoring how well it fits.

**Exhaustiveness** — how hard the search works per compound.

**Grid box / search box** — the 3D region of the receptor that is searched.

**Heavy atoms** — non-hydrogen atoms; ZINC22's measure of molecular size.

**Ligand** — a candidate small molecule.

**logP** — how greasy a molecule is. Most oral drugs sit between about −1 and
5.

**MPI** — the system that lets VinaLC use many cores (and machines) at once.

**PDBQT** — the AutoDock file format: coordinates plus atom types and partial
charges.

**Pose** — one predicted orientation of a ligand in the site.

**Receptor** — the target protein.

**SDF** — a common chemical file format, usually 2D; must be converted before
docking.

**Tranche** — one batch of compounds, processed and checkpointed as a unit.

**VinaLC** — the parallel AutoDock Vina used for docking.

---

## 12. How it works inside

For whoever maintains this next.

| File | Responsibility |
|---|---|
| `pipeline.py` | The five stages, tranche status, stop/resume, ranking, pose export. Importable and runnable as a CLI. |
| `io_parse.py` | Reads and validates `setup.txt`, the tranche list, receptors and boxes. All format knowledge lives here. |
| `zinc22.py` | ZINC22 tranche picking and download-script generation. |
| `chembl.py` | ChEMBL filtering, weight binning and download-script generation. |
| `ligprep.py` | SDF → 3D PDBQT conversion via Open Babel. |
| `resources.py` | Machine detection and the CPU/memory plan. |
| `gui.py` | `PipelineAPI`: every GUI action, as plain functions returning JSON-able dicts. |
| `gui_tk.py` | The Tk window — a thin client of `PipelineAPI` over HTTP. |
| `web_server.py` | Serves `PipelineAPI` on localhost with a token. |
| `gui.sh` | Starts the detached server and opens a window onto it. |

Design rules worth preserving:

- **The GUI never reimplements pipeline logic.** It calls the same functions
  the CLI does. If behaviour must change, change `pipeline.py`.
- **State lives on disk, not in memory.** A tranche's `status.txt` is the
  source of truth, which is what makes runs resumable and lets the window
  recover results after a restart.
- **A stage that is interrupted records nothing**, so it reruns cleanly.
- **The run outlives the window.** The server owns the work.
- **Never call a slow endpoint from the GUI thread.** Anything touching the
  network or a subprocess goes through `DockingFlowApp.call_async`, whose
  worker thread touches no Tk widget and hands its result back through a
  queue that the GUI thread drains. A synchronous call freezes the entire
  window until the far end answers, which with a catalogue outage means
  minutes. A test enforces this.

Run the tests after any change:

```bash
python3 -m unittest discover -s tests -v
```

They cover the whole pipeline offline, using a stand-in docking binary and a
stand-in Open Babel, so they need neither MPI nor a real docking install.

---

## 13. Known limitations

Honest list, as of this writing:

- **Stopping during docking loses that tranche's docking progress.** VinaLC
  cannot resume partway through a tranche, so that tranche restarts. Finished
  tranches are never redone.
- **The browser fallback window is behind the Tk window.** It lacks the ligand
  presets, grid-box buttons, docking settings and pose export.
- **ChEMBL screens depend on Open Babel's behaviour**, some of which is not
  documented by Open Babel itself (how it delimits molecules in a
  multi-molecule PDBQT, and whether it preserves compound names). DockingFlow
  handles both possibilities and verifies compound identity itself, and
  **Check Open Babel** reports what your install actually does — run it before
  trusting a ChEMBL screen.
- **Scores are rankings, not measurements.** See
  [Reading the scores](#reading-the-scores).
- **No automatic disk-space check** before a large screen; use
  `server_check.sh`.
