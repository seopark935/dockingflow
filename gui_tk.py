#!/usr/bin/env python3
"""DockingFlow GUI as a lightweight Tkinter window, for use over X11 (e.g. MobaXterm).

A browser window (Firefox, Chromium, or pywebview's embedded one) forwarded
over X11 is slow: it ships the rendered page to your screen as pixels. Tk
sends small drawing commands instead, so it stays responsive over an SSH
connection. Tkinter ships with Python (on some Linux installs as the
`python3-tk` package).

This window is only a client of the GUI server (`gui.py --web`, started by
`gui.sh`): every button calls the same `PipelineAPI` method the web page
does, over HTTP on localhost. The run itself lives in the server process, so
closing this window never stops it.

Usage (normally via `bash gui.sh`):
    python3 gui_tk.py "http://localhost:8765/?token=..."
"""
from __future__ import annotations

import json
import sys
import tkinter as tk
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from tkinter import filedialog, messagebox, ttk
from typing import Any, Dict, List, Optional

POLL_MS = 1000
MAX_HIT_ROWS = 500  # the full list is in top_hits_combined.tsv

COLORS = {"ok": "#15803d", "err": "#b91c1c", "info": "#1d4ed8", "": "#374151"}

LARGE_SCREEN = 5_000_000  # molecules; above this, confirm before building the list

CHARGES = {"J": "-4", "K": "-3", "L": "-2", "M": "-1", "N": "0", "O": "+1", "P": "+2", "Q": "+3", "R": "+4"}

# (setup.txt key, label, hint) for the Docking tab's settings fields.
SETTINGS = [
    ("filter_percent", "Keep top %", "share of each tranche's ranked ligands kept as hits"),
    ("exhaustiveness", "Exhaustiveness", "search effort per ligand (default 8)"),
    ("num_modes", "Poses per ligand", "binding poses reported (default 9)"),
    ("energy_range", "Energy range", "kcal/mol window of reported poses (default 3)"),
    ("seed", "Random seed", "optional; set for reproducible runs"),
    ("granularity", "Grid spacing (Å)", "optional; default 0.375"),
]


class ApiClient:
    """Calls `PipelineAPI` methods on the GUI server: `api.validate(...)` -> POST /api/validate."""

    def __init__(self, url: str) -> None:
        parts = urllib.parse.urlparse(url)
        self.base = f"{parts.scheme}://{parts.netloc}"
        self.token = (urllib.parse.parse_qs(parts.query).get("token") or [""])[0]

    def __getattr__(self, name: str):
        return lambda *args: self._call(name, list(args))

    def _call(self, name: str, args: list) -> Any:
        req = urllib.request.Request(
            f"{self.base}/api/{name}",
            data=json.dumps(args).encode(),
            method="POST",
            headers={"Content-Type": "application/json", "X-DockingFlow-Token": self.token},
        )
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            with e:
                return json.load(e)
        except (urllib.error.URLError, OSError) as e:
            raise ConnectionError(str(e))


class DockingFlowApp(tk.Tk):
    def __init__(self, api: ApiClient) -> None:
        super().__init__()
        self.api = api
        self.title("DockingFlow")
        self.geometry("1200x800")
        self.minsize(900, 600)
        ttk.Style(self).theme_use("clam")

        self.v = {name: tk.StringVar(self) for name in ("setup_path", "map_path", "workdir", "vinalc_bin", "mpirun_bin")}
        self.no_mpirun = tk.BooleanVar(self, value=False)

        # Ligands tab: ZINC22 tranche picker
        # Defaults: a small drug-like slice, so a first click can't start a billion-molecule screen.
        self.hac_min, self.hac_max = tk.IntVar(self, value=17), tk.IntVar(self, value=17)
        self.logp_min, self.logp_max = tk.DoubleVar(self, value=1.0), tk.DoubleVar(self, value=1.5)
        self.charges = {c: tk.BooleanVar(self, value=c == "N") for c in CHARGES}
        self._zinc_job = None

        # Docking tab: settings saved into setup.txt
        self.settings = {k: tk.StringVar(self) for k, _, _ in SETTINGS}
        self.cores = tk.IntVar(self, value=1)
        self.memory_gb = tk.IntVar(self, value=1)

        self.machine: Optional[Dict[str, Any]] = None
        self.recommended: Optional[Dict[str, Any]] = None
        self.targets: List[Dict[str, Any]] = []  # rows of StringVars
        self.targets_dirty = False
        self._plan_job = None
        self._poll_job = None
        self._last_status: Dict[str, Any] = {}

        self._build()
        self.after(50, self._startup)

    # ---------- API plumbing ----------

    def call(self, name: str, *args: Any) -> Any:
        """Call the server; on a dropped connection, say so and return None."""
        try:
            return getattr(self.api, name)(*args)
        except ConnectionError:
            self.set_message(
                "Lost connection to the DockingFlow GUI server. Any run keeps going on the server; "
                "run 'bash gui.sh' again to reconnect.", "err",
            )
            return None

    def set_message(self, text: str, kind: str = "") -> None:
        self.message.configure(text=text, foreground=COLORS.get(kind, COLORS[""]))

    # ---------- layout ----------

    def _build(self) -> None:
        panes = ttk.PanedWindow(self, orient=tk.HORIZONTAL)
        panes.pack(fill=tk.BOTH, expand=True, padx=8, pady=8)

        left = ttk.Frame(panes, width=440)
        right = ttk.Frame(panes)
        panes.add(left, weight=0)
        panes.add(right, weight=1)

        tabs = ttk.Notebook(left)
        tabs.pack(fill=tk.BOTH, expand=True)
        tabs.add(self._build_ligands(tabs), text="Ligands")
        tabs.add(self._build_targets_tab(tabs), text="Targets")
        tabs.add(self._build_docking(tabs), text="Docking")
        tabs.add(self._build_resources(tabs), text="Resources")

        actions = ttk.Frame(left, padding=(0, 8, 0, 0))
        actions.pack(fill=tk.X)
        ttk.Button(actions, text="Validate", command=self.validate).pack(side=tk.LEFT)
        self.run_btn = ttk.Button(actions, text="Run pipeline", command=self.start_run)
        self.run_btn.pack(side=tk.LEFT, padx=6)
        ttk.Button(actions, text="Nuke workdir", command=self.nuke).pack(side=tk.RIGHT)
        ttk.Button(actions, text="Clean outputs", command=self.clean).pack(side=tk.RIGHT, padx=6)

        self.message = ttk.Label(left, text="Loading...", wraplength=420, justify=tk.LEFT, padding=(2, 8))
        self.message.pack(fill=tk.X)

        self._build_status(right)

    def _file_row(self, parent: ttk.Frame, row: int, label: str, var: tk.StringVar, kind: Optional[str]) -> None:
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=(6, 0), columnspan=2)
        ttk.Entry(parent, textvariable=var).grid(row=row + 1, column=0, sticky="ew")
        if kind:
            ttk.Button(parent, text="Browse", command=lambda: self.browse(var, kind)).grid(
                row=row + 1, column=1, padx=(4, 0))

    def _build_ligands(self, parent: ttk.Notebook) -> ttk.Frame:
        f = ttk.Frame(parent, padding=10)
        f.columnconfigure(0, weight=1)
        pick = ttk.LabelFrame(f, text="Choose ZINC22 tranches", padding=8)
        pick.grid(row=0, column=0, columnspan=2, sticky="ew")

        def range_row(row: int, label: str, lo: tk.Variable, hi: tk.Variable, frm: float, to: float, inc: float) -> None:
            ttk.Label(pick, text=label).grid(row=row, column=0, sticky="w", pady=2)
            for col, var in ((1, lo), (3, hi)):
                sb = ttk.Spinbox(pick, from_=frm, to=to, increment=inc, textvariable=var, width=6,
                                 command=self.schedule_zinc_summary)
                sb.grid(row=row, column=col, padx=2)
                sb.bind("<KeyRelease>", lambda _e: self.schedule_zinc_summary())
            ttk.Label(pick, text="to").grid(row=row, column=2)

        range_row(0, "Heavy atoms", self.hac_min, self.hac_max, 4, 29, 1)
        range_row(1, "logP", self.logp_min, self.logp_max, -5, 9, 0.1)
        ttk.Label(pick, text="Charges").grid(row=2, column=0, sticky="nw", pady=(4, 0))
        charges = ttk.Frame(pick)
        charges.grid(row=2, column=1, columnspan=4, sticky="w", pady=(4, 0))
        for k, (c, label) in enumerate(CHARGES.items()):
            ttk.Checkbutton(charges, text=label, variable=self.charges[c], command=self.schedule_zinc_summary).grid(
                row=k // 5, column=k % 5, sticky="w", padx=(0, 6))
        self.zinc_label = ttk.Label(pick, text="", wraplength=380, justify=tk.LEFT, padding=(0, 6, 0, 0))
        self.zinc_label.grid(row=3, column=0, columnspan=5, sticky="w")
        ttk.Button(pick, text="Create tranche list", command=self.zinc_build).grid(
            row=4, column=0, columnspan=2, sticky="w", pady=(6, 0))

        self._file_row(f, 1, "Tranche list (filled in by the button above)", self.v["map_path"], "file")
        ttk.Button(f, text="Or import a CartBlanche22 download file...", command=self.import_zinc).grid(
            row=3, column=0, sticky="w", pady=(4, 0))
        self._file_row(f, 4, "Work directory (downloads and results; needs lots of disk)", self.v["workdir"], "folder")
        return f

    def _build_docking(self, parent: ttk.Notebook) -> ttk.Frame:
        f = ttk.Frame(parent, padding=10)
        f.columnconfigure(0, weight=1)
        self._file_row(f, 0, "Docking binary (vinalc, or its full path)", self.v["vinalc_bin"], "file")
        self._file_row(f, 2, "MPI launcher", self.v["mpirun_bin"], None)
        ttk.Checkbutton(f, text="Skip MPI launcher (test mock only)", variable=self.no_mpirun).grid(
            row=4, column=0, sticky="w", pady=(6, 0))

        box = ttk.LabelFrame(f, text="Docking settings (saved to the settings file on Validate/Run)", padding=8)
        box.grid(row=5, column=0, columnspan=2, sticky="ew", pady=(10, 0))
        for r, (key, label, hint) in enumerate(SETTINGS):
            ttk.Label(box, text=label).grid(row=r, column=0, sticky="w", pady=2)
            e = ttk.Entry(box, textvariable=self.settings[key], width=8, justify=tk.RIGHT)
            e.grid(row=r, column=1, padx=6)
            ttk.Label(box, text=hint, foreground=COLORS[""]).grid(row=r, column=2, sticky="w")
        self.settings["exhaustiveness"].trace_add("write", lambda *_: self.schedule_plan())

        self._file_row(f, 6, "Settings file (setup.txt)", self.v["setup_path"], "file")
        self.v["setup_path"].trace_add("write", lambda *_: self.after_idle(self._setup_changed))
        return f

    def _build_targets_tab(self, parent: ttk.Notebook) -> ttk.Frame:
        f = ttk.Frame(parent, padding=10)
        ttk.Label(
            f, wraplength=400, justify=tk.LEFT,
            text="Each receptor (a prepared .pdbqt) is docked within its grid box: the box center and its "
                 "size along x/y/z, in Ångströms. Keep boxes around the binding pocket (≤ ~30 Å).",
        ).pack(fill=tk.X)
        self.targets_frame = ttk.Frame(f)
        self.targets_frame.pack(fill=tk.X, pady=6)
        buttons = ttk.Frame(f)
        buttons.pack(fill=tk.X)
        ttk.Button(buttons, text="+ Add receptor", command=self.add_target).pack(side=tk.LEFT)
        ttk.Button(buttons, text="Save targets", command=self.save_targets).pack(side=tk.LEFT, padx=6)
        return f

    def _build_resources(self, parent: ttk.Notebook) -> ttk.Frame:
        f = ttk.Frame(parent, padding=10)
        buttons = ttk.Frame(f)
        buttons.pack(fill=tk.X)
        ttk.Button(buttons, text="Detect this machine", command=self.detect_resources).pack(side=tk.LEFT)
        ttk.Button(buttons, text="Load server report...", command=self.load_report).pack(side=tk.LEFT, padx=6)

        self.specs = ttk.Label(f, text="Detecting...", justify=tk.LEFT, padding=(0, 8))
        self.specs.pack(fill=tk.X)

        ttk.Label(f, text="CPU cores to use").pack(anchor="w")
        self.cores_scale = tk.Scale(f, from_=1, to=1, orient=tk.HORIZONTAL, variable=self.cores,
                                    command=lambda _: self.schedule_plan(), highlightthickness=0)
        self.cores_scale.pack(fill=tk.X)
        ttk.Label(f, text="Memory budget (GB)").pack(anchor="w", pady=(6, 0))
        self.mem_scale = tk.Scale(f, from_=1, to=1, orient=tk.HORIZONTAL, variable=self.memory_gb,
                                  command=lambda _: self.schedule_plan(), highlightthickness=0)
        self.mem_scale.pack(fill=tk.X)

        self.plan_label = ttk.Label(f, text="", justify=tk.LEFT, wraplength=400, padding=(0, 8))
        self.plan_label.pack(fill=tk.X)

        buttons2 = ttk.Frame(f)
        buttons2.pack(fill=tk.X)
        ttk.Button(buttons2, text="Use recommended", command=self.use_recommended).pack(side=tk.LEFT)
        return f

    def _tree(self, parent: tk.Widget, columns: List[tuple], height: int) -> ttk.Treeview:
        frame = ttk.Frame(parent)
        frame.pack(fill=tk.BOTH, expand=True)
        tree = ttk.Treeview(frame, columns=[c[0] for c in columns], show="headings", height=height)
        for name, width in columns:
            tree.heading(name, text=name)
            tree.column(name, width=width, anchor="w")
        scroll = ttk.Scrollbar(frame, orient=tk.VERTICAL, command=tree.yview)
        tree.configure(yscrollcommand=scroll.set)
        tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)
        return tree

    def _build_status(self, right: ttk.Frame) -> None:
        box = ttk.LabelFrame(right, text="Tranches", padding=6)
        box.pack(fill=tk.X)
        self.tranche_tree = self._tree(
            box, [("Label", 110), ("LogP", 60), ("Size", 80), ("Status", 130), ("Error", 380)], 5)
        self.tranche_tree.tag_configure("DONE", foreground=COLORS["ok"])
        self.tranche_tree.tag_configure("FAILED", foreground=COLORS["err"])

        box = ttk.LabelFrame(right, text="Log", padding=6)
        box.pack(fill=tk.BOTH, pady=8)
        self.log = tk.Text(box, height=10, wrap=tk.WORD, state=tk.DISABLED, font=("TkFixedFont", 9))
        scroll = ttk.Scrollbar(box, orient=tk.VERTICAL, command=self.log.yview)
        self.log.configure(yscrollcommand=scroll.set)
        self.log.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)

        box = ttk.LabelFrame(right, text="Top hits (combined, ranked by affinity)", padding=6)
        box.pack(fill=tk.BOTH, expand=True)
        self.hits_tree = self._tree(
            box, [("Ligand", 220), ("Affinity (kcal/mol)", 130), ("Receptor", 140), ("Tranche", 100)], 10)

    # ---------- startup ----------

    def _startup(self) -> None:
        d = self.call("get_defaults")
        if d is None:
            return
        for key in self.v:
            self.v[key].set(d[key])
        self.detect_resources()
        self.schedule_zinc_summary()
        s = self.call("get_status")
        if s:
            if s["phase"] == "running":
                self.set_message("A run is in progress on the server; showing its live status.", "info")
            elif s["phase"] in ("done", "error"):
                self.set_message(s["message"], "ok" if s["phase"] == "done" else "err")
            else:
                self.set_message("Work through the tabs left to right, then click Validate.")
        self.poll()

    def _setup_changed(self) -> None:
        self.load_targets()
        self.load_settings()
        self.schedule_plan()

    # ---------- settings <-> setup.txt ----------

    def load_settings(self) -> None:
        if not self.v["setup_path"].get():
            return
        res = self.call("get_settings", self.v["setup_path"].get())
        if not res or not res["ok"]:
            return
        saved = res["settings"]
        for key, _, _ in SETTINGS:
            self.settings[key].set(saved.get(key, ""))
        self._saved_budget = (saved.get("cores"), saved.get("memory_gb"))
        self.apply_saved_budget()

    def apply_saved_budget(self) -> None:
        """Show the settings file's saved CPU/memory budget, if it fits this machine."""
        cores, mem = getattr(self, "_saved_budget", (None, None))
        if not self.machine or not cores:
            return
        try:
            if int(cores) <= self.machine["logical_cpus"]:
                self.cores.set(int(cores))
            if mem and self.machine.get("mem_total_gb") and float(mem) <= self.machine["mem_total_gb"]:
                self.memory_gb.set(int(float(mem)))
        except ValueError:
            return
        self.schedule_plan()

    def save_settings(self) -> bool:
        """Write every docking setting and the CPU/memory budget into setup.txt."""
        values = {key: self.settings[key].get() for key, _, _ in SETTINGS}
        cores, mem = self.budget()
        if cores is not None:
            values["cores"] = str(cores)
            values["memory_gb"] = str(int(mem)) if mem else ""
        res = self.call("save_settings", self.v["setup_path"].get(), values)
        if res is None:
            return False
        if not res["ok"]:
            self.set_message(f"Docking settings: {res['message']}", "err")
        return res["ok"]

    def save_all(self) -> bool:
        """Save targets (if edited) and settings; False if anything was invalid."""
        if self.targets_dirty and not self.save_targets():
            return False
        return self.save_settings()

    # ---------- ZINC22 tranche picker ----------

    def _zinc_filters(self) -> Optional[Dict[str, Any]]:
        try:
            return {
                "hac_min": int(self.hac_min.get()), "hac_max": int(self.hac_max.get()),
                "logp_min": float(self.logp_min.get()), "logp_max": float(self.logp_max.get()),
                "charges": [c for c, v in self.charges.items() if v.get()],
            }
        except (tk.TclError, ValueError):
            return None

    def schedule_zinc_summary(self) -> None:
        if self._zinc_job:
            self.after_cancel(self._zinc_job)
        self._zinc_job = self.after(400, self.refresh_zinc_summary)

    def refresh_zinc_summary(self) -> None:
        self._zinc_job = None
        filters = self._zinc_filters()
        if filters is None:
            self.zinc_label.configure(text="Enter numbers for the ranges.")
            return
        self.zinc_label.configure(text="Counting (the first time downloads ZINC22's index, ~10 MB)...")
        self.update_idletasks()
        res = self.call("zinc_summary", filters)
        if res is None:
            return
        if not res["ok"]:
            self.zinc_label.configure(text=res["message"])
            return
        self.zinc_label.configure(
            text=f"{res['tranches']:,} tranche(s), about {res['molecules']:,} molecules."
            + ("" if res["molecules"] < LARGE_SCREEN else "  That's a very large screen; consider narrowing it.")
        )

    def zinc_build(self) -> None:
        filters = self._zinc_filters()
        if filters is None:
            self.set_message("Enter numbers for the heavy atom and logP ranges.", "err")
            return
        summary = self.call("zinc_summary", filters)
        if summary and summary.get("ok") and summary["molecules"] > LARGE_SCREEN and not messagebox.askyesno(
            "Large screen",
            f"This selection is about {summary['molecules']:,} molecules in {summary['tranches']:,} tranches.\n\n"
            "Docking that many takes a very long time and a lot of disk. Create the list anyway?",
        ):
            return
        self.set_message("Getting the file list from ZINC22...", "info")
        self.update_idletasks()
        res = self.call("zinc_build", filters)
        if res is None:
            return
        if res["ok"]:
            self.v["map_path"].set(res["map_path"])
        self.set_message(res["message"], "ok" if res["ok"] else "err")

    # ---------- file pickers + ZINC import ----------

    def browse(self, var: tk.StringVar, kind: str) -> None:
        start = Path(var.get() or self.v["setup_path"].get()).expanduser()
        initial = str(start if start.is_dir() else start.parent)
        path = (filedialog.askdirectory(initialdir=initial) if kind == "folder"
                else filedialog.askopenfilename(initialdir=initial))
        if path:
            var.set(path)

    def import_zinc(self) -> None:
        initial = str(Path(self.v["map_path"].get() or ".").expanduser().parent)
        path = filedialog.askopenfilename(
            initialdir=initial, title="ZINC22 downloader file from CartBlanche22",
            filetypes=[("Downloader scripts", "*.curl *.wget *.txt"), ("All files", "*")])
        if not path:
            return
        res = self.call("import_zinc_downloader", path)
        if res is None:
            return
        if res["ok"]:
            self.v["map_path"].set(res["map_path"])
        self.set_message(res["message"], "ok" if res["ok"] else "err")

    # ---------- docking targets ----------

    def load_targets(self) -> None:
        if not self.v["setup_path"].get():
            return
        res = self.call("get_targets", self.v["setup_path"].get())
        if res is None:
            return
        if not res["ok"]:
            self.set_message(f"Couldn't read docking targets: {res['message']}", "err")
            return
        self.targets = [self._target_vars(t) for t in res["targets"]] or [self._target_vars(None)]
        self.targets_dirty = False
        self.render_targets()

    def _target_vars(self, t: Optional[Dict[str, Any]]) -> Dict[str, Any]:
        t = t or {"receptor": "", "center": ["", "", ""], "size": ["20", "20", "20"]}
        row = {
            "receptor": tk.StringVar(self, value=t["receptor"]),
            "center": [tk.StringVar(self, value=v) for v in t["center"]],
            "size": [tk.StringVar(self, value=v) for v in t["size"]],
        }
        for var in [row["receptor"]] + row["center"] + row["size"]:
            var.trace_add("write", lambda *_: setattr(self, "targets_dirty", True))
        return row

    def render_targets(self) -> None:
        for child in self.targets_frame.winfo_children():
            child.destroy()
        for i, row in enumerate(self.targets):
            box = ttk.LabelFrame(self.targets_frame, text=f"Receptor {i + 1}", padding=6)
            box.pack(fill=tk.X, pady=3)
            box.columnconfigure(1, weight=1)
            ttk.Entry(box, textvariable=row["receptor"]).grid(row=0, column=0, columnspan=4, sticky="ew")
            ttk.Button(box, text="Browse", command=lambda v=row["receptor"]: self.browse(v, "file")).grid(
                row=0, column=4, padx=(4, 0))
            ttk.Button(box, text="Remove", command=lambda k=i: self.remove_target(k)).grid(row=0, column=5, padx=(4, 0))
            for r, (label, key) in enumerate((("Center", "center"), ("Size Å", "size")), start=1):
                ttk.Label(box, text=label, width=7).grid(row=r, column=0, sticky="w", pady=(4, 0))
                for k, var in enumerate(row[key]):
                    ttk.Entry(box, textvariable=var, width=10, justify=tk.RIGHT).grid(
                        row=r, column=1 + k, sticky="ew", padx=2, pady=(4, 0))

    def add_target(self) -> None:
        self.targets.append(self._target_vars(None))
        self.targets_dirty = True
        self.render_targets()

    def remove_target(self, i: int) -> None:
        del self.targets[i]
        self.targets_dirty = True
        self.render_targets()

    def save_targets(self) -> bool:
        payload = [{
            "receptor": row["receptor"].get().strip(),
            "center": [v.get().strip() for v in row["center"]],
            "size": [v.get().strip() for v in row["size"]],
        } for row in self.targets]
        res = self.call("save_targets", self.v["setup_path"].get(), payload)
        if res is None:
            return False
        self.set_message(res["message"], "err" if not res["ok"] else ("info" if res["warnings"] else "ok"))
        if res["ok"]:
            self.targets_dirty = False
            self.schedule_plan()  # box size changes memory per worker
        return res["ok"]

    # ---------- resources ----------

    def detect_resources(self) -> None:
        self.show_machine(self.call("detect_resources"))

    def load_report(self) -> None:
        path = filedialog.askopenfilename(title="server_check.sh output (server_report.txt)")
        if path:
            self.show_machine(self.call("load_server_report", path))

    def show_machine(self, res: Optional[Dict[str, Any]]) -> None:
        if not res:
            return
        if not res["ok"]:
            self.set_message(res["message"], "err")
            return
        self.machine, self.recommended = res["resources"], res["recommended"]
        m = self.machine
        lines = [f"{m['hostname']}  ({m['source']})",
                 f"CPU: {m['physical_cores']} physical cores, {m['logical_cpus']} logical"
                 + (f", load {m['load1']:.1f}" if m.get("load1") is not None else "")]
        if m.get("mem_total_gb"):
            lines.append(f"RAM: {m['mem_total_gb']:.0f} GB total"
                         + (f", {m['mem_available_gb']:.0f} GB available" if m.get("mem_available_gb") else ""))
        lines += self.recommended["notes"]
        self.specs.configure(text="\n".join(lines))
        self.cores_scale.configure(to=m["logical_cpus"])
        self.mem_scale.configure(to=max(1, int(m.get("mem_total_gb") or 1)),
                                 state=tk.NORMAL if m.get("mem_total_gb") else tk.DISABLED)
        self.use_recommended()
        self.apply_saved_budget()

    def use_recommended(self) -> None:
        if not self.recommended:
            return
        self.cores.set(self.recommended["cores"])
        if self.recommended.get("memory_gb"):
            self.memory_gb.set(int(self.recommended["memory_gb"]))
        self.schedule_plan()

    def budget(self) -> tuple:
        if not self.machine:
            return None, None
        return self.cores.get(), (float(self.memory_gb.get()) if self.machine.get("mem_total_gb") else None)

    def schedule_plan(self) -> None:
        if self._plan_job:
            self.after_cancel(self._plan_job)
        self._plan_job = self.after(250, self.refresh_plan)

    def refresh_plan(self) -> None:
        self._plan_job = None
        cores, mem = self.budget()
        if cores is None or not self.v["setup_path"].get():
            return
        try:
            exh = int(self.settings["exhaustiveness"].get() or 0) or None
        except ValueError:
            exh = None
        p = self.call("plan_budget", self.v["setup_path"].get(), cores, mem, exh)
        if p is None:
            return
        if not p["ok"]:
            self.plan_label.configure(text=f"Can't plan yet: {p['message']}")
            return
        lines = [
            f"VinaLC: mpirun -np {p['mpi_ranks']}  (1 master + {p['workers']} workers x ~{p['exhaustiveness']} threads)",
            f"Est. memory: {p['est_memory_gb']:.1f} GB  (~{p['rank_memory_gb']:.1f} GB per worker for your grid box)",
            f"Limited by: {p['limited_by']}",
        ]
        if cores > self.machine["physical_cores"]:
            lines.append("More cores than physical cores: hyperthreads add little for Vina.")
        if self.machine.get("mem_available_gb") and p["est_memory_gb"] > self.machine["mem_available_gb"]:
            lines.append("Estimate exceeds currently available RAM.")
        self.plan_label.configure(text="\n".join(lines))

    # ---------- validate / run / clean ----------

    def validate(self) -> None:
        if not self.save_all():
            return
        self.set_message("Validating...", "info")
        self.update_idletasks()
        res = self.call("validate", self.v["setup_path"].get(), self.v["map_path"].get(), self.v["workdir"].get())
        if res:
            self.set_message(res["message"], "ok" if res["ok"] else "err")

    def start_run(self) -> None:
        if not self.save_all():
            return
        res = self.call(
            "start_run", self.v["setup_path"].get(), self.v["map_path"].get(), self.v["workdir"].get(),
            self.v["vinalc_bin"].get() or "vinalc", self.v["mpirun_bin"].get() or "mpirun",
            self.no_mpirun.get(),
        )
        if res is None:
            return
        if not res["ok"]:
            self.set_message(res["message"], "err")
            return
        self.set_message("Run started. You can close this window; the run continues on the server.", "info")
        self.poll()

    def clean(self) -> None:
        wd = self.v["workdir"].get()
        if messagebox.askyesno("Clean outputs", f"Delete generated tranche outputs under:\n{wd}\n\nThe workdir itself is kept."):
            res = self.call("clean", wd)
            if res:
                self.set_message("Cleaned tranche outputs." if res["ok"] else res["message"], "ok" if res["ok"] else "err")

    def nuke(self) -> None:
        wd = self.v["workdir"].get()
        if messagebox.askyesno("Nuke workdir", f"DELETE THE ENTIRE WORKDIR:\n{wd}\n\nThis cannot be undone. Continue?",
                               icon=messagebox.WARNING):
            res = self.call("nuke", wd)
            if res:
                self.set_message("Workdir deleted." if res["ok"] else res["message"], "ok" if res["ok"] else "err")

    # ---------- live status ----------

    def poll(self) -> None:
        if self._poll_job:
            self.after_cancel(self._poll_job)
            self._poll_job = None
        s = self.call("get_status")
        if s is not None:
            self.render_status(s)
        self._poll_job = self.after(POLL_MS, self.poll)

    def render_status(self, s: Dict[str, Any]) -> None:
        last = self._last_status
        if s["tranches"] != last.get("tranches"):
            self.tranche_tree.delete(*self.tranche_tree.get_children())
            for t in s["tranches"]:
                tag = "DONE" if t["status"] == "DONE" else ("FAILED" if t["status"].startswith("FAILED") else "")
                self.tranche_tree.insert("", tk.END, values=(t["label"], t["log_p"], t["size"], t["status"], t["error"] or ""),
                                         tags=(tag,))

        if s["log"] != last.get("log"):
            at_bottom = self.log.yview()[1] > 0.98
            self.log.configure(state=tk.NORMAL)
            self.log.delete("1.0", tk.END)
            self.log.insert(tk.END, "\n".join(s["log"]))
            self.log.configure(state=tk.DISABLED)
            if at_bottom:
                self.log.see(tk.END)

        if s["top_hits"] != last.get("top_hits"):
            self.hits_tree.delete(*self.hits_tree.get_children())
            for h in s["top_hits"][:MAX_HIT_ROWS]:
                self.hits_tree.insert("", tk.END, values=(h["ligand"], f"{h['affinity']:.3f}", h["receptor"], h["tranche"]))

        if s["phase"] != last.get("phase"):
            self.run_btn.configure(state=tk.DISABLED if s["phase"] == "running" else tk.NORMAL)
            if last and s["phase"] == "done":
                self.set_message(s["message"], "ok")
            elif last and s["phase"] == "error":
                self.set_message(s["message"], "err")
        self._last_status = s


def main() -> int:
    if len(sys.argv) != 2 or "token=" not in sys.argv[1]:
        print(__doc__)
        return 2
    app = DockingFlowApp(ApiClient(sys.argv[1]))
    app.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
