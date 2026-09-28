#!/usr/bin/env python3
"""DockingFlow desktop GUI.

This module wraps the existing, unmodified pipeline (`pipeline.py` /
`io_parse.py`) in a small desktop application so a run can be configured and
watched without touching a terminal.

Architecture
------------
`pywebview` hosts a single HTML/CSS/JS page (`gui_assets/index.html`) inside a
native OS window. Python exposes a `PipelineAPI` instance as the window's
`js_api`; the frontend calls its methods as `pywebview.api.<method>(...)`,
which pywebview marshals to/from JSON automatically (Python dicts/lists/
primitives <-> JS objects/arrays/values).

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

import threading
import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import webview

import io_parse
import pipeline
import resources

REPO_ROOT = Path(__file__).resolve().parent


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
        self._window: webview.Window | None = None

    def set_window(self, window: webview.Window) -> None:
        """Wired up once, right after the window is created, for file dialogs."""
        self._window = window

    # ---- defaults / file pickers ----

    def get_defaults(self) -> dict[str, str]:
        """Prefill the form with this repo's own config files, for convenience."""
        return {
            "setup_path": str(REPO_ROOT / "setup.txt"),
            "map_path": str(REPO_ROOT / "tranches.txt"),
            "workdir": str(REPO_ROOT / "run"),
            "vinalc_bin": "vinalc",
            "mpirun_bin": "mpirun",
        }

    def browse_file(self) -> str | None:
        if not self._window:
            return None
        result = self._window.create_file_dialog(webview.OPEN_DIALOG)
        return result[0] if result else None

    def browse_folder(self) -> str | None:
        if not self._window:
            return None
        result = self._window.create_file_dialog(webview.FOLDER_DIALOG)
        return result[0] if result else None

    # ---- machine resources / CPU + memory budget ----

    def detect_resources(self) -> dict[str, Any]:
        """Specs of the machine the GUI (and so the pipeline) is running on, plus a suggested budget."""
        try:
            res = resources.detect_local()
            return {"ok": True, "resources": res, "recommended": resources.recommend(res)}
        except Exception as exc:
            return {"ok": False, "message": str(exc)}

    def load_server_report(self) -> dict[str, Any] | None:
        """Pick a `server_check.sh` output file and size the run for that machine instead."""
        path = self.browse_file()
        if not path:
            return None
        try:
            res = resources.parse_server_report(Path(path).read_text(errors="replace"))
            res["source"] = f"server report ({Path(path).name})"
            return {"ok": True, "resources": res, "recommended": resources.recommend(res)}
        except Exception as exc:
            return {"ok": False, "message": str(exc)}

    def plan_budget(self, setup_path: str, cores: int, memory_gb: float | None) -> dict[str, Any]:
        """What a core/memory budget means for VinaLC (ranks, per-rank memory), for live display."""
        try:
            return {"ok": True, **resources.plan(setup_path, int(cores), memory_gb)}
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
    api = PipelineAPI()
    html_path = REPO_ROOT / "gui_assets" / "index.html"
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
