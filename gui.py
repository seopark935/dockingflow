#!/usr/bin/env python3
"""DockingFlow GUI: a desktop window, or a web page served from the docking server.

This module wraps the existing pipeline (`pipeline.py` / `io_parse.py`) in a
small app so a run can be configured and watched without touching a
terminal. It runs in one of two modes, sharing the same frontend
(`gui_assets/index.html`) and the same `PipelineAPI`:

  - Desktop (`python3 gui.py`): `pywebview` hosts the page in a native OS
    window and exposes `PipelineAPI` as the window's `js_api`; the frontend
    calls `pywebview.api.<method>(...)`, which pywebview marshals to/from
    JSON.
  - Web (`python3 gui.py --web`): for a headless server. `web_server.py`
    serves the page and a JSON endpoint per `PipelineAPI` method on
    localhost, and a small shim in the page routes the same
    `pywebview.api.<method>(...)` calls over HTTP. View it from a laptop
    through an SSH tunnel. Standard library only — pywebview isn't needed.

Architecture
------------
The actual docking pipeline runs on a background thread (`PipelineAPI._run`)
so the window stays responsive while downloads/docking are in progress. The
frontend polls `get_status()` every ~700ms and re-renders; there is no push
channel from Python to JS, by design, to keep the threading model simple
(one writer thread, stateless reads from the polling loop).

A run isolates per-tranche failures: if one tranche's download or docking
fails, its error is recorded and the loop continues to the next tranche
rather than aborting the whole run. This differs slightly from the CLI's
`pipeline.py main()`, which stops the whole run on the first failure — the
GUI favors showing partial progress over an all-or-nothing run.
"""
from __future__ import annotations

import argparse
import math
import os
import shutil
import threading
import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import io_parse
import pipeline
import resources
import zinc22
import zinc_split

REPO_ROOT = Path(__file__).resolve().parent

# The GUI's own working copy of the config files. Git ignores it, so editing
# settings in the GUI never blocks a `git pull` of code updates. Seeded on
# first launch from the example files checked into the repo.
PROJECT_DIR = REPO_ROOT / "project"
PROJECT_SEED_FILES = ["setup.txt", "recList.txt", "geoList.txt", "protein.pdbqt"]

# Above this box side (Å), docking is effectively "blind" (whole-protein):
# fine, but slower, and it needs more exhaustiveness to search reliably.
LARGE_BOX_ANGSTROM = 30


def ensure_project_dir() -> Path:
    """Create `project/` from the repo's example config on first use; never overwrite it."""
    if not (PROJECT_DIR / "setup.txt").exists():
        PROJECT_DIR.mkdir(exist_ok=True)
        for name in PROJECT_SEED_FILES:
            if not (PROJECT_DIR / name).exists():
                shutil.copy2(REPO_ROOT / name, PROJECT_DIR / name)
        # No CPU/memory budget yet: the GUI recommends one for this machine,
        # and saves the user's choice on the first Validate/Run. Docking
        # settings start at VinaLC's standard defaults.
        resources.set_setup_keys(str(PROJECT_DIR / "setup.txt"), {
            "cores": None, "memory_gb": None, "energy_range": "3", "filter_percent": "10",
        })
    return PROJECT_DIR


@dataclass
class TrancheProgress:
    """One row of the frontend's tranche table."""

    label: str
    log_p: float
    size: str
    status: str = "PENDING"
    error: str | None = None


@dataclass
class RunState:
    """Snapshot of an in-progress or finished run, as seen by the frontend.

    A fresh instance replaces the old one at the start of every run, so
    stale state from a previous run never leaks into a new one.
    """

    phase: str = "idle"  # idle | running | done | error
    message: str = ""
    log: list[str] = field(default_factory=list)
    tranches: list[TrancheProgress] = field(default_factory=list)
    top_hits: list[dict[str, Any]] = field(default_factory=list)

    def log_line(self, text: str) -> None:
        self.log.append(text)
        if len(self.log) > 500:
            del self.log[:-500]  # bound memory on a long-running screen


class PipelineAPI:
    """Methods exposed to the JS frontend as `pywebview.api.<name>(...)`.

    Every public method here is called directly from JavaScript, so each one
    returns a JSON-serializable value (dict/list/str/bool/number/None) and
    never raises — failures are caught and reported back as
    `{"ok": False, "message": "..."}` so the frontend can display them
    without pywebview's own error dialog.
    """

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._state = RunState()
        self._window: Any = None  # webview.Window, desktop mode only

    def set_window(self, window: Any) -> None:
        """Wired up once, right after the window is created, for file dialogs."""
        self._window = window

    # ---- defaults / file pickers ----

    def get_defaults(self) -> dict[str, str]:
        """Prefill the form with the GUI's `project/` config (created on first launch)."""
        project = ensure_project_dir()
        return {
            "setup_path": str(project / "setup.txt"),
            "map_path": str(project / "zinc22_tranches.txt"),  # written by the Ligands tab
            "workdir": str(project / "run"),
            "vinalc_bin": "vinalc",
            "mpirun_bin": "mpirun",
        }

    def browse_file(self) -> str | None:
        """Native file dialog (desktop mode). The web frontend uses `list_dir` instead."""
        if not self._window:
            return None
        import webview

        result = self._window.create_file_dialog(webview.OPEN_DIALOG)
        return result[0] if result else None

    def browse_folder(self) -> str | None:
        """Native folder dialog (desktop mode). The web frontend uses `list_dir` instead."""
        if not self._window:
            return None
        import webview

        result = self._window.create_file_dialog(webview.FOLDER_DIALOG)
        return result[0] if result else None

    def list_dir(self, path: str) -> dict[str, Any]:
        """Directory listing for the web frontend's in-page file picker.

        `path` may be a file (lists its directory) or empty (lists the repo).
        """
        try:
            p = Path(path).expanduser() if path else REPO_ROOT
            if not p.is_absolute():
                p = REPO_ROOT / p
            p = p.resolve()
            while not p.is_dir():  # a file, or a path that doesn't exist yet
                p = p.parent
            entries = []
            with os.scandir(p) as it:
                for e in it:
                    try:
                        entries.append({"name": e.name, "is_dir": e.is_dir()})
                    except OSError:
                        continue
            entries.sort(key=lambda e: (not e["is_dir"], e["name"].startswith("."), e["name"].lower()))
            return {"ok": True, "path": str(p), "parent": str(p.parent), "entries": entries}
        except Exception as exc:
            return {"ok": False, "message": str(exc)}

    # ---- ZINC import ----

    def import_zinc_downloader(self, downloader_path: str) -> dict[str, Any]:
        """Split a CartBlanche22 downloader file into per-tranche scripts + a tranche map.

        Writes `project/zinc22_scripts/<code>.curl` and
        `project/zinc22_tranches.txt` (the same place as the Ligands tab's
        picker); the frontend then points the "Tranches map" field at the
        latter.
        """
        try:
            src = Path(downloader_path).expanduser().resolve()
            project = ensure_project_dir()
            out_dir = project / "zinc22_scripts"
            map_path = project / "zinc22_tranches.txt"
            counts = zinc_split.split_to_files(src, out_dir, map_path)
            n_files = sum(counts.values())
            return {
                "ok": True,
                "map_path": str(map_path),
                "message": (
                    f"Imported {src.name}: {len(counts)} tranche(s), {n_files} download(s).\n"
                    f"Tranches: {', '.join(counts)}\nTranche map: {map_path}"
                ),
            }
        except Exception as exc:
            return {"ok": False, "message": str(exc)}

    # ---- ZINC22 tranche picker: choices -> tranche map, no CartBlanche website ----

    @staticmethod
    def _zinc_selection(filters: dict[str, Any]) -> list[dict[str, Any]]:
        index = zinc22.load_index(ensure_project_dir() / "zinc22_index.json")
        return zinc22.select(
            index,
            int(filters["hac_min"]), int(filters["hac_max"]),
            float(filters["logp_min"]), float(filters["logp_max"]),
            list(filters["charges"]),
        )

    def zinc_summary(self, filters: dict[str, Any]) -> dict[str, Any]:
        """How much of ZINC22 the current choices cover (first call downloads the index, ~10 MB)."""
        try:
            return {"ok": True, **zinc22.summarize(self._zinc_selection(filters))}
        except Exception as exc:
            return {"ok": False, "message": f"Couldn't read the ZINC22 index: {exc}"}

    def zinc_build(self, filters: dict[str, Any]) -> dict[str, Any]:
        """Build the per-tranche download scripts + tranche map for the chosen ZINC22 tranches."""
        try:
            project = ensure_project_dir()
            map_path = project / "zinc22_tranches.txt"
            res = zinc22.build_tranche_map(self._zinc_selection(filters), project / "zinc22_scripts", map_path)
            return {
                "ok": True,
                "map_path": str(map_path),
                "message": f"Tranche list ready: {res['tranches']} tranche(s), {res['files']} file(s) to download.",
            }
        except Exception as exc:
            return {"ok": False, "message": str(exc)}

    # ---- docking settings <-> setup.txt ----

    SETTING_KEYS = ("filter_percent", "exhaustiveness", "num_modes", "energy_range", "seed", "granularity",
                    "cores", "memory_gb")

    def get_settings(self, setup_path: str) -> dict[str, Any]:
        """setup.txt's docking settings, as strings for the form ("" = not set / default)."""
        try:
            setup = io_parse.load_setup(Path(setup_path))
            return {"ok": True, "settings": {k: setup.get(k, "") for k in self.SETTING_KEYS}}
        except Exception as exc:
            return {"ok": False, "message": str(exc)}

    def save_settings(self, setup_path: str, settings: dict[str, Any]) -> dict[str, Any]:
        """Validate the form's docking settings and write them into setup.txt.

        Blank optional fields are removed from setup.txt so the defaults
        apply. Nothing is written unless every value is valid.
        """
        try:
            updates = {k: (str(settings[k]).strip() or None) for k in self.SETTING_KEYS if k in settings}
            for required in ("filter_percent", "cores"):
                if required in updates and updates[required] is None:
                    raise ValueError(f"{required} can't be blank")

            candidate = io_parse.load_setup(Path(setup_path))
            for k, v in updates.items():
                if v is None:
                    candidate.pop(k, None)
                else:
                    candidate[k] = v
            try:
                fp = float(candidate.get("filter_percent", "nan"))
            except ValueError:
                raise ValueError("Keep top % must be a number")
            if not 0 < fp <= 100:
                raise ValueError("Keep top % must be between 0 and 100")
            # Same checks the pipeline itself runs (skipping the memory-based rank
            # planning, which needs valid docking targets).
            io_parse.load_vinalc_options({k: v for k, v in candidate.items() if k != "memory_gb"})
            if candidate.get("memory_gb"):
                io_parse._positive_number(candidate, "memory_gb")

            resources.set_setup_keys(setup_path, updates)
            return {"ok": True, "message": f"Saved settings to {Path(setup_path).name}."}
        except Exception as exc:
            return {"ok": False, "message": str(exc)}

    # ---- docking targets (recList.txt + geoList.txt) ----

    @staticmethod
    def _target_list_paths(setup_path: str) -> tuple[Path, Path]:
        setup = io_parse.load_setup(Path(setup_path))
        setup_dir = Path(setup_path).expanduser().resolve().parent
        rec = Path(setup["recList"]) if "recList" in setup else setup_dir / "recList.txt"
        geo = Path(setup["geoList"]) if "geoList" in setup else setup_dir / "geoList.txt"
        return rec, geo

    def get_targets(self, setup_path: str) -> dict[str, Any]:
        """The receptor + grid box rows from the setup's recList/geoList, for editing.

        Lenient on purpose (unlike `io_parse.load_docking_targets`): missing
        receptors or malformed lines are shown so they can be fixed here.
        """
        try:
            rec, geo = self._target_list_paths(setup_path)
            rec_lines = [ln.strip() for ln in rec.read_text().splitlines() if ln.strip()] if rec.exists() else []
            geo_lines = [ln.split() for ln in geo.read_text().splitlines() if ln.strip()] if geo.exists() else []
            targets = []
            for i in range(max(len(rec_lines), len(geo_lines))):
                receptor = str((rec.parent / rec_lines[i]).resolve()) if i < len(rec_lines) else ""
                box = geo_lines[i] if i < len(geo_lines) else []
                box = (box + [""] * 6)[:6]
                targets.append({"receptor": receptor, "center": box[:3], "size": box[3:]})
            return {"ok": True, "targets": targets, "rec_list": str(rec), "geo_list": str(geo)}
        except Exception as exc:
            return {"ok": False, "message": str(exc)}

    def fit_box_to_receptor(self, receptor: str) -> dict[str, Any]:
        """Blind docking: a box around the whole receptor (5 Å margin), each side capped at 80 Å."""
        try:
            box = io_parse.fit_box(io_parse.read_coordinates(Path(receptor).expanduser()), padding=5)
            msg = (f"Box fitted to the whole receptor ({box['atoms']:,} atoms, "
                   f"{' x '.join(f'{e:g}' for e in box['extent'])} Å across) for blind docking.")
            if box["capped"]:
                msg += (f" The protein is wider than {io_parse.MAX_BOX_SIDE:g} Å along {', '.join(box['capped'])}, "
                        "so the box was capped there and misses some of its surface. If you know the binding "
                        "site, center the box on it instead.")
            return {"ok": True, **box, "message": msg}
        except Exception as exc:
            return {"ok": False, "message": str(exc)}

    def fit_box_to_site(self, reference: str) -> dict[str, Any]:
        """Known site: a box around a reference ligand (e.g. a co-crystallized inhibitor), >= 20 Å per side."""
        try:
            box = io_parse.fit_box(io_parse.read_coordinates(Path(reference).expanduser()), padding=5, min_side=20)
            return {"ok": True, **box, "message": (
                f"Box centered on {Path(reference).name} ({box['atoms']} atoms) with a 5 Å margin "
                f"(at least 20 Å per side). Make sure that file is in the same coordinate frame as your receptor, "
                f"e.g. both taken from the same PDB entry.")}
        except Exception as exc:
            return {"ok": False, "message": str(exc)}

    def save_targets(self, setup_path: str, targets: list[dict[str, Any]]) -> dict[str, Any]:
        """Validate the edited receptor/grid-box rows and write recList.txt + geoList.txt.

        Each target is `{"receptor": path, "center": [x, y, z], "size": [x, y, z]}`.
        Nothing is written unless every row is valid.
        """
        try:
            if not targets:
                raise ValueError("Add at least one receptor.")
            rec_rows, geo_rows, warnings = [], [], []
            for i, t in enumerate(targets, start=1):
                receptor = Path(str(t.get("receptor", "")).strip()).expanduser()
                if not str(receptor) or str(receptor) == ".":
                    raise ValueError(f"Target {i}: choose a receptor file.")
                receptor = receptor.resolve()
                if not receptor.is_file():
                    raise ValueError(f"Target {i}: receptor not found: {receptor}")
                if receptor.suffix.lower() != ".pdbqt":
                    raise ValueError(f"Target {i}: receptor must be a prepared .pdbqt file: {receptor.name}")
                problem, placeholder = io_parse.check_receptor(receptor)
                if problem:
                    raise ValueError(f"Target {i}: {receptor.name}: {problem}")
                if placeholder:
                    warnings.append(f"Target {i}: {receptor.name} is the placeholder receptor; "
                                    "choose your real prepared receptor before a real screen.")

                try:
                    center = [float(v) for v in t.get("center", [])]
                    size = [float(v) for v in t.get("size", [])]
                except (TypeError, ValueError):
                    raise ValueError(f"Target {i}: box center and size must be numbers")
                if len(center) != 3 or len(size) != 3 or not all(map(math.isfinite, center + size)):
                    raise ValueError(f"Target {i}: box needs center x/y/z and size x/y/z")
                if min(size) <= 0:
                    raise ValueError(f"Target {i}: box sizes must be > 0")
                if max(size) > io_parse.MAX_BOX_SIDE:
                    raise ValueError(
                        f"Target {i}: box sides can be at most {io_parse.MAX_BOX_SIDE:g} Å (blind docking of a "
                        f"typical protein). Use 'Fit to whole protein' to size it automatically."
                    )
                if max(size) > LARGE_BOX_ANGSTROM:
                    warnings.append(
                        f"Target {i}: a {max(size):g} Å box is blind docking (the whole protein). That's fine, "
                        f"but slower; consider exhaustiveness 16-32 in the Docking tab so the larger space is "
                        f"searched thoroughly."
                    )
                rec_rows.append(str(receptor))
                geo_rows.append(" ".join(f"{v:g}" for v in center + size))

            rec, geo = self._target_list_paths(setup_path)
            rec.write_text("\n".join(rec_rows) + "\n")
            geo.write_text("\n".join(geo_rows) + "\n")
            resources.set_setup_keys(setup_path, {"recList": str(rec), "geoList": str(geo)})
            message = f"Saved {len(rec_rows)} docking target(s) to {rec.name} / {geo.name}."
            return {"ok": True, "message": "\n".join([message] + warnings), "warnings": warnings}
        except Exception as exc:
            return {"ok": False, "message": str(exc)}

    # ---- machine resources / CPU + memory budget ----

    def detect_resources(self) -> dict[str, Any]:
        """Specs of the machine the GUI (and so the pipeline) is running on, plus a suggested budget."""
        try:
            res = resources.detect_local()
            return {"ok": True, "resources": res, "recommended": resources.recommend(res)}
        except Exception as exc:
            return {"ok": False, "message": str(exc)}

    def load_server_report(self, path: str | None = None) -> dict[str, Any] | None:
        """Size the run from a `server_check.sh` output file instead of this machine.

        `path` comes from the web frontend's picker; in desktop mode it's
        omitted and a native file dialog is shown.
        """
        path = path or self.browse_file()
        if not path:
            return None
        try:
            res = resources.parse_server_report(Path(path).read_text(errors="replace"))
            res["source"] = f"server report ({Path(path).name})"
            return {"ok": True, "resources": res, "recommended": resources.recommend(res)}
        except Exception as exc:
            return {"ok": False, "message": str(exc)}

    def plan_budget(
        self, setup_path: str, cores: int, memory_gb: float | None, exhaustiveness: int | None = None
    ) -> dict[str, Any]:
        """What a core/memory budget means for VinaLC (ranks, per-rank memory), for live display.

        `exhaustiveness` lets the form preview an unsaved value.
        """
        try:
            return {"ok": True, **resources.plan(setup_path, int(cores), memory_gb, exhaustiveness)}
        except Exception as exc:
            return {"ok": False, "message": str(exc)}

    def save_budget(self, setup_path: str, cores: int, memory_gb: float | None) -> dict[str, Any]:
        """Write the chosen budget into setup.txt (`cores=`, `memory_gb=`), e.g. for a CLI run on the server."""
        try:
            resources.write_budget_to_setup(setup_path, int(cores), memory_gb)
            return {"ok": True, "message": f"Saved cores={int(cores)}, memory_gb={memory_gb} to {setup_path}"}
        except Exception as exc:
            return {"ok": False, "message": str(exc)}

    @staticmethod
    def _load_setup_with_budget(setup_path: str, cores: int | None, memory_gb: float | None) -> dict[str, str]:
        """setup.txt, with the GUI's core/memory budget (if one is set) taking precedence."""
        setup = io_parse.load_setup(Path(setup_path))
        if cores:
            setup["cores"] = str(int(cores))
        if memory_gb:
            setup["memory_gb"] = str(memory_gb)
        return setup

    # ---- validation (read-only, safe to call anytime) ----

    def validate(
        self,
        setup_path: str,
        map_path: str,
        workdir: str,
        cores: int | None = None,
        memory_gb: float | None = None,
    ) -> dict[str, Any]:
        """Run the same checks `pipeline.py` runs before Stage 0, without starting anything."""
        try:
            workdir_p = Path(workdir).expanduser().resolve()
            workdir_p.mkdir(parents=True, exist_ok=True)
            setup = self._load_setup_with_budget(setup_path, cores, memory_gb)
            tranches = io_parse.load_tranches_tsv(Path(map_path))
            io_parse.validate_inputs(setup, tranches, workdir_p)
            return {"ok": True, "message": io_parse.format_summary(setup, tranches, workdir_p)}
        except Exception as exc:
            return {"ok": False, "message": str(exc)}

    # ---- run orchestration ----

    def start_run(
        self,
        setup_path: str,
        map_path: str,
        workdir: str,
        vinalc_bin: str,
        mpirun_bin: str,
        no_mpirun: bool,
        cores: int | None = None,
        memory_gb: float | None = None,
    ) -> dict[str, Any]:
        """Kick off a full pipeline run on a background thread; returns immediately.

        `cores`/`memory_gb`, when given (from the Resources panel), override
        setup.txt's values for this run only.
        """
        with self._lock:
            if self._state.phase == "running":
                return {"ok": False, "message": "A run is already in progress."}
            self._state = RunState(phase="running", message="Starting...")

        thread = threading.Thread(
            target=self._run,
            args=(setup_path, map_path, workdir, vinalc_bin, mpirun_bin, no_mpirun, cores, memory_gb),
            daemon=True,
        )
        thread.start()
        return {"ok": True}

    def _run(
        self,
        setup_path: str,
        map_path: str,
        workdir: str,
        vinalc_bin: str,
        mpirun_bin: str,
        no_mpirun: bool,
        cores: int | None,
        memory_gb: float | None,
    ) -> None:
        """The actual pipeline run, executed on a background thread.

        Mutates `self._state` in place as it progresses; `get_status()` reads
        that same object from the main thread. There's no explicit lock
        around individual field writes here — CPython's GIL makes each
        attribute assignment atomic enough for a status display that's
        re-rendered from scratch every poll, and `start_run` already
        prevents two runs from writing to the same `RunState` concurrently.
        """
        state = self._state
        try:
            workdir_p = Path(workdir).expanduser().resolve()
            workdir_p.mkdir(parents=True, exist_ok=True)

            setup = self._load_setup_with_budget(setup_path, cores, memory_gb)
            tranches = io_parse.load_tranches_tsv(Path(map_path))
            io_parse.validate_inputs(setup, tranches, workdir_p)
            targets = io_parse.load_docking_targets(setup)
            options = io_parse.load_vinalc_options(setup)

            state.tranches = [
                TrancheProgress(pipeline.tranche_label(t), t.log_p, io_parse.tranche_size(t)) for t in tranches
            ]
            state.log_line(f"Loaded {len(tranches)} tranche(s).")
            if "mpi_ranks" in setup:
                state.log_line(f"Note: setup.txt sets mpi_ranks={setup['mpi_ranks']}, which overrides the CPU/memory budget.")
            state.log_line(
                f"Budget: cores={setup['cores']}, memory_gb={setup.get('memory_gb', 'unlimited')} "
                f"-> mpirun -np {options.mpi_ranks}"
            )

            (workdir_p / "tranches").mkdir(exist_ok=True)
            tdirs = []
            for t, prog in zip(tranches, state.tranches):
                tdir = pipeline.materialize_tranche(workdir_p, t, Path(setup_path))
                tdirs.append(tdir)
                prog.status = pipeline.read_status(tdir)

            filter_percent = float(setup["filter_percent"])
            mpirun = None if no_mpirun else mpirun_bin

            for tdir, prog in zip(tdirs, state.tranches):
                try:
                    state.log_line(f"[{prog.label}] downloading...")
                    pipeline.download_tranche(tdir)
                    prog.status = pipeline.read_status(tdir)

                    state.log_line(f"[{prog.label}] unpacking...")
                    pipeline.unpack_tranche(tdir)
                    prog.status = pipeline.read_status(tdir)

                    state.log_line(f"[{prog.label}] docking...")
                    pipeline.dock_tranche(
                        tdir,
                        targets,
                        options=options,
                        vinalc_bin=vinalc_bin,
                        mpirun_bin=mpirun,
                    )
                    prog.status = pipeline.read_status(tdir)

                    state.log_line(f"[{prog.label}] ranking...")
                    pipeline.rank_and_filter_tranche(tdir, filter_percent)
                    prog.status = pipeline.read_status(tdir)
                    state.log_line(f"[{prog.label}] done.")
                except Exception as exc:
                    prog.status = pipeline.read_status(tdir)
                    prog.error = str(exc)
                    state.log_line(f"[{prog.label}] FAILED: {exc}")
                    continue  # isolate failures: one bad tranche shouldn't kill the whole run

            combined_path = pipeline.combine_results(workdir_p, tdirs)
            state.top_hits = self._read_top_hits(combined_path)
            state.phase = "done"
            state.message = f"Finished. Combined results: {combined_path}"
        except Exception as exc:
            state.phase = "error"
            state.message = str(exc)
            state.log_line(f"ERROR: {exc}\n{traceback.format_exc()}")

    def tranche_report(self, workdir: str, label: str) -> dict[str, Any]:
        """Everything useful for debugging one tranche, as one readable text report.

        Status and counts, VinaLC's command and any errors it printed, and the
        ends of the download/docking logs — the same files listed in README's
        Troubleshooting section, gathered in one place.
        """
        try:
            tdir = (Path(workdir).expanduser().resolve() / "tranches" / label)
            if not tdir.is_dir():
                raise ValueError(f"No such tranche folder: {tdir}")
            logs, dock = tdir / "logs", tdir / "docking"

            def read(path: Path) -> str:
                return path.read_text(errors="replace").strip() if path.exists() else ""

            def count_lines(path: Path) -> int:
                return max(0, sum(1 for _ in open(path, errors="replace")) - 1) if path.exists() else 0

            def section(title: str, body: str) -> str:
                return f"==== {title} ====\n{body.strip() or '(empty)'}\n"

            n_archives = pipeline.count_ligand_archives(tdir / "download") if (tdir / "download").exists() else 0
            counts = [
                f"status:              {pipeline.read_status(tdir)}",
                f"downloaded archives: {n_archives}",
                f"failed downloads:    {len(read(tdir / 'download_failures.txt').splitlines())}",
                f"molecules unpacked:  {count_lines(tdir / 'ligand_index.tsv')}",
                read(dock / "dock_summary.txt").replace("jobs=", "docking jobs:        ").replace("results=", "docked with a pose:  "),
                f"unique ligands ranked: {count_lines(tdir / 'results' / 'ranked_all.tsv')} (ZINC repeats some compounds across archives)",
            ]
            vina_log = dock / f"{pipeline.VINALC_REC_LIST}_{pipeline.VINALC_LIG_LIST}.log.gz"
            vina_log_head = ""
            if vina_log.exists():
                import gzip

                with gzip.open(vina_log, "rt", errors="replace") as f:
                    vina_log_head = "".join(line for _, line in zip(range(40), f))

            parts = [
                f"Tranche {label}  ({tdir})\n",
                section("Summary", "\n".join(c for c in counts if c.strip())),
                section("VinaLC errors", "\n".join(pipeline.vinalc_errors(logs / "dock.stderr.log", logs / "dock.stdout.log", limit=10))),
                section("Docking command (docking/command.txt)", read(dock / "command.txt") + "\n" + read(dock / "dock_timing.txt")),
                section("Docking stderr, last 40 lines (logs/dock.stderr.log)", "\n".join(pipeline.tail_lines(logs / "dock.stderr.log"))),
                section("Docking stdout, last 15 lines (logs/dock.stdout.log)", "\n".join(pipeline.tail_lines(logs / "dock.stdout.log", 15))),
                section("VinaLC per-ligand log, first 40 lines (docking/*.log.gz)", vina_log_head),
                section("Download failures (download_failures.txt)", read(tdir / "download_failures.txt")),
                section("Download stderr, last 15 lines (logs/download.stderr.log)", "\n".join(pipeline.tail_lines(logs / "download.stderr.log", 15))),
                section("Unpack warnings (unpack_warnings.txt)", read(tdir / "unpack_warnings.txt")),
            ]
            return {"ok": True, "report": "\n".join(parts), "folder": str(tdir)}
        except Exception as exc:
            return {"ok": False, "message": str(exc)}

    @staticmethod
    def _read_top_hits(path: Path) -> list[dict[str, Any]]:
        if not path.exists():
            return []
        rows = []
        for line in path.read_text().splitlines()[1:]:
            if not line.strip():
                continue
            lig, aff, receptor, tranche = line.split("\t")
            rows.append({"ligand": lig, "affinity": float(aff), "receptor": receptor, "tranche": tranche})
        return rows

    # ---- polling ----

    def get_status(self) -> dict[str, Any]:
        """Called every ~700ms by the frontend to re-render the current run state."""
        state = self._state
        return {
            "phase": state.phase,
            "message": state.message,
            "log": state.log[-200:],
            "tranches": [
                {
                    "label": t.label,
                    "log_p": t.log_p,
                    "size": t.size,
                    "status": t.status,
                    "error": t.error,
                }
                for t in state.tranches
            ],
            "top_hits": state.top_hits,
        }

    # ---- destructive actions (confirmed client-side before being called) ----

    def clean(self, workdir: str) -> dict[str, Any]:
        try:
            pipeline.clean_workdir(Path(workdir).expanduser())
            return {"ok": True}
        except Exception as exc:
            return {"ok": False, "message": str(exc)}

    def nuke(self, workdir: str) -> dict[str, Any]:
        try:
            pipeline.nuke_workdir(Path(workdir).expanduser())
            return {"ok": True}
        except Exception as exc:
            return {"ok": False, "message": str(exc)}


def main() -> None:
    parser = argparse.ArgumentParser(description="DockingFlow GUI")
    parser.add_argument(
        "--web", action="store_true",
        help="Serve the GUI as a web page (for a headless server) instead of opening a desktop window",
    )
    parser.add_argument("--port", type=int, default=8765, help="Web mode: port to listen on (default 8765)")
    parser.add_argument(
        "--host", default="127.0.0.1",
        help="Web mode: address to bind (default 127.0.0.1, reachable only through an SSH tunnel)",
    )
    parser.add_argument(
        "--viewer", metavar="URL",
        help="Open a desktop window showing an already-running web GUI (used by gui.sh over X11)",
    )
    args = parser.parse_args()

    if args.viewer:
        # Just a window onto the web server: closing it never stops a run,
        # because the run lives in the server process. Remote X11 displays
        # (e.g. MobaXterm) rarely support GPU rendering, so turn it off.
        os.environ.setdefault("QTWEBENGINE_CHROMIUM_FLAGS", "--disable-gpu")
        os.environ.setdefault("WEBKIT_DISABLE_COMPOSITING_MODE", "1")
        import webview

        webview.create_window("DockingFlow", args.viewer, width=1200, height=820, min_size=(860, 600),
                              background_color="#0f1115")
        webview.start()
        return

    api = PipelineAPI()
    html_path = REPO_ROOT / "gui_assets" / "index.html"

    if args.web:
        import web_server

        web_server.serve(api, html_path, host=args.host, port=args.port)
        return

    import webview

    window = webview.create_window(
        "DockingFlow",
        str(html_path),
        js_api=api,
        width=1100,
        height=780,
        min_size=(860, 600),
        background_color="#0f1115",
    )
    api.set_window(window)
    webview.start()


if __name__ == "__main__":
    main()
