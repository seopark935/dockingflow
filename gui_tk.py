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
import tkinter.font as tkfont
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from tkinter import filedialog, messagebox, ttk
from typing import Any, Dict, List, Optional

POLL_MS = 1000
MAX_HIT_ROWS = 500  # the full list is in top_hits_combined.tsv

# Light theme. Plain colors only (no images), so it renders the same over X11.
P = {
    "bg": "#eef0f3", "card": "#ffffff", "border": "#d5d9e0", "text": "#111827", "muted": "#6b7280",
    "button": "#e8eaee", "button_hover": "#dde0e6",
    "accent": "#2563eb", "accent_hover": "#1d4ed8", "accent_disabled": "#9dbaf3",
    "ok": "#15803d", "ok_bg": "#dcfce7", "err": "#b91c1c", "err_bg": "#fee2e2",
    "info": "#1d4ed8", "info_bg": "#dbeafe", "neutral_bg": "#e5e7eb",
    "stripe": "#f7f8fa", "select": "#cfe0fd", "log_bg": "#111827", "log_fg": "#d1d5db",
}
COLORS = {"ok": P["ok"], "err": P["err"], "info": P["info"], "": P["muted"]}
MESSAGE_BG = {"ok": P["ok_bg"], "err": P["err_bg"], "info": P["info_bg"], "": P["neutral_bg"]}

# Friendlier names for the pipeline's tranche statuses.
STATUS_TEXT = {
    "PENDING": "Waiting", "INIT": "Queued", "DOWNLOADED": "Downloaded", "UNPACKED": "Unpacked",
    "DOCKED": "Docked", "DONE": "Done",
}
PHASE_PILL = {  # run phase -> (text, fg, bg)
    "idle": ("Idle", P["muted"], P["neutral_bg"]), "running": ("Running", P["info"], P["info_bg"]),
    "done": ("Finished", P["ok"], P["ok_bg"]), "error": ("Error", P["err"], P["err_bg"]),
}

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
        self.geometry("1280x840")
        self.minsize(980, 640)
        self._style()

        self.v = {name: tk.StringVar(self) for name in ("setup_path", "map_path", "workdir", "vinalc_bin", "mpirun_bin")}
        self.no_mpirun = tk.BooleanVar(self, value=False)

        # Ligands tab: ZINC22 tranche picker
        # Default: one ~40k-molecule tranche, a sensible pilot run; widen from there.
        self.hac_min, self.hac_max = tk.IntVar(self, value=12), tk.IntVar(self, value=12)
        self.logp_min, self.logp_max = tk.DoubleVar(self, value=1.0), tk.DoubleVar(self, value=1.0)
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
        self.message.configure(text=text, fg=P["text"] if not kind else COLORS[kind], bg=MESSAGE_BG.get(kind, P["neutral_bg"]))

    # ---------- look ----------

    def _style(self) -> None:
        base = tkfont.nametofont("TkDefaultFont")
        base.configure(size=10)
        for name in ("TkTextFont", "TkMenuFont", "TkHeadingFont"):
            tkfont.nametofont(name).configure(family=base.cget("family"), size=10)
        self.font_bold = base.copy()
        self.font_bold.configure(weight="bold")
        self.font_title = base.copy()
        self.font_title.configure(size=15, weight="bold")
        self.font_mono = tkfont.nametofont("TkFixedFont").copy()
        self.font_mono.configure(size=9)

        self.configure(background=P["bg"])
        st = ttk.Style(self)
        st.theme_use("clam")
        # Content widgets sit on white cards by default; the window frame is gray ("App.*").
        st.configure(".", background=P["card"], foreground=P["text"], bordercolor=P["border"],
                     lightcolor=P["card"], darkcolor=P["card"], troughcolor=P["bg"], focuscolor=P["accent"])
        st.configure("App.TFrame", background=P["bg"])
        st.configure("App.TLabel", background=P["bg"])
        st.configure("Muted.TLabel", foreground=P["muted"])
        st.configure("AppMuted.TLabel", background=P["bg"], foreground=P["muted"])
        st.configure("Title.TLabel", background=P["bg"], font=self.font_title)
        st.configure("Section.TLabel", font=self.font_bold)

        st.configure("TLabelframe", background=P["card"], bordercolor=P["border"], relief="solid", borderwidth=1)
        st.configure("TLabelframe.Label", background=P["card"], foreground=P["muted"], font=self.font_bold)
        st.configure("Card.TFrame", background=P["card"], bordercolor=P["border"], relief="solid", borderwidth=1)

        st.configure("TNotebook", background=P["bg"], borderwidth=0, tabmargins=(0, 0, 0, 0))
        st.configure("TNotebook.Tab", background=P["button"], foreground=P["muted"], padding=(16, 7),
                     bordercolor=P["border"])
        st.map("TNotebook.Tab", background=[("selected", P["card"])], foreground=[("selected", P["text"])],
               expand=[("selected", (0, 0, 0, 0))])

        st.configure("TButton", background=P["button"], padding=(12, 5), bordercolor=P["border"], relief="flat")
        st.map("TButton", background=[("active", P["button_hover"]), ("disabled", P["button"])],
               foreground=[("disabled", P["muted"])])
        st.configure("Accent.TButton", background=P["accent"], foreground="#ffffff", bordercolor=P["accent"],
                     font=self.font_bold)
        st.map("Accent.TButton", background=[("active", P["accent_hover"]), ("disabled", P["accent_disabled"])],
               foreground=[("disabled", "#ffffff")])
        st.configure("Danger.TButton", foreground=P["err"])
        st.map("Danger.TButton", background=[("active", P["err_bg"])])

        for widget in ("TEntry", "TSpinbox"):
            st.configure(widget, fieldbackground="#ffffff", bordercolor=P["border"], padding=4,
                         lightcolor=P["border"], darkcolor=P["border"])
            st.map(widget, bordercolor=[("focus", P["accent"])], lightcolor=[("focus", P["accent"])])
        st.configure("TCheckbutton", background=P["card"])
        st.map("TCheckbutton", background=[("active", P["card"])])

        st.configure("Treeview", background=P["card"], fieldbackground=P["card"], rowheight=24,
                     bordercolor=P["border"], borderwidth=0)
        st.map("Treeview", background=[("selected", P["select"])], foreground=[("selected", P["text"])])
        st.configure("Treeview.Heading", background=P["stripe"], foreground=P["muted"], font=self.font_bold,
                     relief="flat", padding=(6, 4))
        st.map("Treeview.Heading", background=[("active", P["button"])])
        st.configure("Accent.Horizontal.TProgressbar", background=P["accent"], troughcolor=P["neutral_bg"],
                     bordercolor=P["neutral_bg"], lightcolor=P["accent"], darkcolor=P["accent"], thickness=8)
        st.configure("TPanedwindow", background=P["bg"])
        st.configure("Vertical.TScrollbar", background=P["button"], troughcolor=P["card"], bordercolor=P["card"],
                     arrowcolor=P["muted"])

    # ---------- layout ----------

    def _build(self) -> None:
        outer = ttk.Frame(self, style="App.TFrame", padding=(14, 10, 14, 12))
        outer.pack(fill=tk.BOTH, expand=True)

        header = ttk.Frame(outer, style="App.TFrame")
        header.pack(fill=tk.X, pady=(0, 10))
        ttk.Label(header, text="DockingFlow", style="Title.TLabel").pack(side=tk.LEFT)
        self.subtitle = ttk.Label(header, text="virtual screening with VinaLC", style="AppMuted.TLabel")
        self.subtitle.pack(side=tk.LEFT, padx=(10, 0), pady=(5, 0))
        self.pill = tk.Label(header, text="Idle", font=self.font_bold, padx=12, pady=3)
        self.pill.pack(side=tk.RIGHT)
        self.set_phase("idle")

        panes = ttk.PanedWindow(outer, orient=tk.HORIZONTAL)
        panes.pack(fill=tk.BOTH, expand=True)
        left = ttk.Frame(panes, style="App.TFrame", width=470)
        right = ttk.Frame(panes, style="App.TFrame")
        panes.add(left, weight=0)
        panes.add(right, weight=1)

        tabs = ttk.Notebook(left)
        tabs.pack(fill=tk.BOTH, expand=True, padx=(0, 12))
        tabs.add(self._build_ligands(tabs), text="Ligands")
        tabs.add(self._build_targets_tab(tabs), text="Targets")
        tabs.add(self._build_docking(tabs), text="Docking")
        tabs.add(self._build_resources(tabs), text="Resources")

        actions = ttk.Frame(left, style="Card.TFrame", padding=10)
        actions.pack(fill=tk.X, padx=(0, 12), pady=(10, 0))
        row = ttk.Frame(actions)
        row.pack(fill=tk.X)
        ttk.Button(row, text="Validate", command=self.validate).pack(side=tk.LEFT)
        self.run_btn = ttk.Button(row, text="Run pipeline", style="Accent.TButton", command=self.start_run)
        self.run_btn.pack(side=tk.LEFT, padx=8)
        ttk.Button(row, text="Nuke workdir", style="Danger.TButton", command=self.nuke).pack(side=tk.RIGHT)
        ttk.Button(row, text="Clean outputs", command=self.clean).pack(side=tk.RIGHT, padx=8)
        self.message = tk.Label(actions, text="Loading...", wraplength=420, justify=tk.LEFT, anchor="w",
                                padx=10, pady=8)
        self.message.pack(fill=tk.X, pady=(10, 0))
        self.set_message("Loading...")

        self._build_status(right)

    def set_phase(self, phase: str) -> None:
        text, fg, bg = PHASE_PILL.get(phase, PHASE_PILL["idle"])
        self.pill.configure(text=text, fg=fg, bg=bg)

    def _file_row(self, parent: ttk.Frame, row: int, label: str, var: tk.StringVar, kind: Optional[str]) -> None:
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=(6, 0), columnspan=2)
        ttk.Entry(parent, textvariable=var).grid(row=row + 1, column=0, sticky="ew")
        if kind:
            ttk.Button(parent, text="Browse", command=lambda: self.browse(var, kind)).grid(
                row=row + 1, column=1, padx=(4, 0))

    def _build_ligands(self, parent: ttk.Notebook) -> ttk.Frame:
        f = ttk.Frame(parent, padding=14)
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
        f = ttk.Frame(parent, padding=14)
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
            ttk.Label(box, text=hint, style="Muted.TLabel").grid(row=r, column=2, sticky="w")
        self.settings["exhaustiveness"].trace_add("write", lambda *_: self.schedule_plan())

        self._file_row(f, 6, "Settings file (setup.txt)", self.v["setup_path"], "file")
        self.v["setup_path"].trace_add("write", lambda *_: self.after_idle(self._setup_changed))
        return f

    def _build_targets_tab(self, parent: ttk.Notebook) -> ttk.Frame:
        f = ttk.Frame(parent, padding=14)
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
        f = ttk.Frame(parent, padding=14)
        buttons = ttk.Frame(f)
        buttons.pack(fill=tk.X)
        ttk.Button(buttons, text="Detect this machine", command=self.detect_resources).pack(side=tk.LEFT)
        ttk.Button(buttons, text="Load server report...", command=self.load_report).pack(side=tk.LEFT, padx=6)

        self.specs = ttk.Label(f, text="Detecting...", justify=tk.LEFT, padding=(0, 8))
        self.specs.pack(fill=tk.X)

        ttk.Label(f, text="CPU cores to use").pack(anchor="w")
        self.cores_scale = tk.Scale(f, from_=1, to=1, orient=tk.HORIZONTAL, variable=self.cores,
                                    command=lambda _: self.schedule_plan(), **self._scale_look())
        self.cores_scale.pack(fill=tk.X)
        ttk.Label(f, text="Memory budget (GB)").pack(anchor="w", pady=(6, 0))
        self.mem_scale = tk.Scale(f, from_=1, to=1, orient=tk.HORIZONTAL, variable=self.memory_gb,
                                  command=lambda _: self.schedule_plan(), **self._scale_look())
        self.mem_scale.pack(fill=tk.X)

        self.plan_label = ttk.Label(f, text="", justify=tk.LEFT, wraplength=400, padding=(0, 8))
        self.plan_label.pack(fill=tk.X)

        buttons2 = ttk.Frame(f)
        buttons2.pack(fill=tk.X)
        ttk.Button(buttons2, text="Use recommended", command=self.use_recommended).pack(side=tk.LEFT)
        return f

    @staticmethod
    def _scale_look() -> Dict[str, Any]:
        return dict(highlightthickness=0, bg=P["card"], fg=P["text"], troughcolor=P["neutral_bg"],
                    activebackground=P["accent_hover"], sliderrelief=tk.FLAT, bd=0, sliderlength=22, width=12)

    def _tree(self, parent: tk.Widget, columns: List[tuple], height: int) -> ttk.Treeview:
        frame = ttk.Frame(parent)
        frame.pack(fill=tk.BOTH, expand=True)
        tree = ttk.Treeview(frame, columns=[c[0] for c in columns], show="headings", height=height)
        for name, width, *anchor in columns:
            tree.heading(name, text=name, anchor=anchor[0] if anchor else "w")
            tree.column(name, width=width, anchor=anchor[0] if anchor else "w", stretch=True)
        tree.tag_configure("odd", background=P["stripe"])
        scroll = ttk.Scrollbar(frame, orient=tk.VERTICAL, command=tree.yview)
        tree.configure(yscrollcommand=scroll.set)
        tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)
        return tree

    def _card(self, parent: tk.Widget, title: str, **pack: Any) -> ttk.Frame:
        card = ttk.Frame(parent, style="Card.TFrame", padding=(12, 10))
        card.pack(**pack)
        head = ttk.Frame(card)
        head.pack(fill=tk.X, pady=(0, 6))
        ttk.Label(head, text=title, style="Section.TLabel").pack(side=tk.LEFT)
        card.head = head
        return card

    def _build_status(self, right: ttk.Frame) -> None:
        box = self._card(right, "Tranches", fill=tk.X)
        ttk.Button(box.head, text="Details...", command=self.show_details).pack(side=tk.RIGHT)
        self.progress_text = ttk.Label(box.head, text="", style="Muted.TLabel")
        self.progress_text.pack(side=tk.RIGHT, padx=10)
        self.progress = ttk.Progressbar(box, style="Accent.Horizontal.TProgressbar", mode="determinate")
        self.progress.pack(fill=tk.X, pady=(0, 8))
        self.tranche_tree = self._tree(
            box, [("Tranche", 100), ("LogP", 60, "e"), ("Size", 80), ("Status", 120), ("Error", 420)], 6)
        self.tranche_tree.tag_configure("DONE", foreground=P["ok"])
        self.tranche_tree.tag_configure("FAILED", foreground=P["err"])
        self.tranche_tree.tag_configure("ACTIVE", foreground=P["info"])
        self.tranche_tree.bind("<Double-1>", lambda _e: self.show_details())

        box = self._card(right, "Log", fill=tk.BOTH, pady=10)
        frame = ttk.Frame(box)
        frame.pack(fill=tk.BOTH, expand=True)
        self.log = tk.Text(frame, height=9, wrap=tk.WORD, state=tk.DISABLED, font=self.font_mono,
                           bg=P["log_bg"], fg=P["log_fg"], relief=tk.FLAT, padx=10, pady=8,
                           highlightthickness=0, insertbackground=P["log_fg"])
        scroll = ttk.Scrollbar(frame, orient=tk.VERTICAL, command=self.log.yview)
        self.log.configure(yscrollcommand=scroll.set)
        self.log.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)

        box = self._card(right, "Top hits", fill=tk.BOTH, expand=True)
        ttk.Label(box.head, text="best predicted binders first (more negative = stronger)",
                  style="Muted.TLabel").pack(side=tk.LEFT, padx=10)
        self.hits_tree = self._tree(
            box, [("Ligand", 240), ("Affinity (kcal/mol)", 140, "e"), ("Receptor", 150), ("Tranche", 110)], 10)

    # ---------- per-tranche debug report ----------

    def show_details(self) -> None:
        sel = self.tranche_tree.selection()
        if not sel:
            self.set_message("Select a tranche in the Tranches table first, then click Details.", "info")
            return
        label = self.tranche_tree.item(sel[0])["values"][0]
        win = tk.Toplevel(self)
        win.title(f"DockingFlow - tranche {label}")
        win.geometry("980x680")
        win.configure(background=P["bg"])
        top = ttk.Frame(win, style="App.TFrame", padding=(12, 10, 12, 6))
        top.pack(fill=tk.X)
        ttk.Label(top, text=f"Tranche {label}", style="Title.TLabel").pack(side=tk.LEFT)
        body = ttk.Frame(win, style="Card.TFrame", padding=6)
        body.pack(fill=tk.BOTH, expand=True, padx=12, pady=(0, 12))
        text = tk.Text(body, wrap=tk.NONE, font=self.font_mono, relief=tk.FLAT, padx=10, pady=8,
                       highlightthickness=0, bg=P["card"], fg=P["text"])
        ys = ttk.Scrollbar(body, orient=tk.VERTICAL, command=text.yview)
        xs = ttk.Scrollbar(body, orient=tk.HORIZONTAL, command=text.xview)
        text.configure(yscrollcommand=ys.set, xscrollcommand=xs.set)
        ys.pack(side=tk.RIGHT, fill=tk.Y)
        xs.pack(side=tk.BOTTOM, fill=tk.X)
        text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        text.tag_configure("head", font=self.font_bold, foreground=P["accent"])
        text.tag_configure("err", foreground=P["err"])

        def load() -> None:
            res = self.call("tranche_report", self.v["workdir"].get(), label)
            text.configure(state=tk.NORMAL)
            text.delete("1.0", tk.END)
            report = res["report"] if res and res.get("ok") else (res or {}).get("message", "Couldn't load the report.")
            for line in report.splitlines():
                tag = "head" if line.startswith("====") else ("err" if pipeline_error(line) else "")
                text.insert(tk.END, line + "\n", tag)
            text.configure(state=tk.DISABLED)

        def copy() -> None:
            self.clipboard_clear()
            self.clipboard_append(text.get("1.0", tk.END))
            self.set_message(f"Copied the report for {label}; paste it wherever you need it.", "ok")

        ttk.Button(top, text="Close", command=win.destroy).pack(side=tk.RIGHT)
        ttk.Button(top, text="Copy report", command=copy).pack(side=tk.RIGHT, padx=8)
        ttk.Button(top, text="Refresh", command=load).pack(side=tk.RIGHT)
        load()

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
        self.subtitle.configure(text=f"virtual screening with VinaLC on {m['hostname']}")
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
            selected = [self.tranche_tree.item(i)["values"][0] for i in self.tranche_tree.selection()]
            self.tranche_tree.delete(*self.tranche_tree.get_children())
            for k, t in enumerate(s["tranches"]):
                status = t["status"]
                if status == "DONE":
                    tag, text = "DONE", "Done"
                elif status.startswith("FAILED"):
                    tag, text = "FAILED", "Failed: " + status.split("_", 1)[-1].lower()
                else:
                    tag, text = ("ACTIVE" if status not in ("PENDING", "INIT") else ""), STATUS_TEXT.get(status, status)
                iid = self.tranche_tree.insert(
                    "", tk.END, values=(t["label"], t["log_p"], t["size"], text, t["error"] or ""),
                    tags=tuple(x for x in (tag, "odd" if k % 2 else "") if x))
                if t["label"] in selected:
                    self.tranche_tree.selection_add(iid)
            n = len(s["tranches"])
            done = sum(1 for t in s["tranches"] if t["status"] == "DONE")
            failed = sum(1 for t in s["tranches"] if t["status"].startswith("FAILED"))
            self.progress.configure(maximum=max(n, 1), value=done + failed)
            self.progress_text.configure(
                text=(f"{done} of {n} done" + (f", {failed} failed" if failed else "")) if n else "")

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
            for k, h in enumerate(s["top_hits"][:MAX_HIT_ROWS]):
                self.hits_tree.insert("", tk.END, values=(h["ligand"], f"{h['affinity']:.2f}", h["receptor"], h["tranche"]),
                                      tags=("odd",) if k % 2 else ())

        if s["phase"] != last.get("phase"):
            self.set_phase(s["phase"])
            self.run_btn.configure(state=tk.DISABLED if s["phase"] == "running" else tk.NORMAL)
            if last and s["phase"] == "done":
                failed = [t["label"] for t in s["tranches"] if t["status"].startswith("FAILED")]
                if failed:
                    self.set_message(f"Finished, but {len(failed)} tranche(s) failed ({', '.join(failed[:5])}"
                                     f"{'...' if len(failed) > 5 else ''}). Select one and click Details to see why.",
                                     "err")
                else:
                    self.set_message(s["message"], "ok")
            elif last and s["phase"] == "error":
                self.set_message(s["message"], "err")
        self._last_status = s


def pipeline_error(line: str) -> bool:
    """Whether a report line looks like an error worth highlighting."""
    low = line.lower()
    return any(w in low for w in ("error", "failed (rc", "failed:", "cannot", "can't", "not found"))

def main() -> int:
    if len(sys.argv) != 2 or "token=" not in sys.argv[1]:
        print(__doc__)
        return 2
    app = DockingFlowApp(ApiClient(sys.argv[1]))
    app.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
