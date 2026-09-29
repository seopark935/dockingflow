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

Written for first-time users: every section has a "?" button explaining it
(texts in `HELP`), every setting starts at a sensible default, and the
ligand library and charges can be picked from common presets.

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

from zinc22 import CHARGE_PRESETS, LIGAND_PRESETS, LOGP_BINS, HAC_MAX, HAC_MIN

POLL_MS = 1000
MAX_HIT_ROWS = 500  # the full list is in top_hits_combined.tsv
LARGE_SCREEN = 5_000_000  # molecules; above this, confirm before building the list
DA_PER_HEAVY_ATOM = 14  # rough average for drug-like molecules, for the "~ Da" hint

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

CHARGES = {"J": "-4", "K": "-3", "L": "-2", "M": "-1", "N": "0", "O": "+1", "P": "+2", "Q": "+3", "R": "+4"}
CUSTOM = "Custom"
DEFAULT_PRESET = "Quick test (~1,000 molecules)"

# (setup.txt key, label, default shown when setup.txt doesn't set it) for the Docking tab.
SETTINGS = [
    ("filter_percent", "Keep top %", "10"),
    ("exhaustiveness", "Exhaustiveness", "8"),
    ("num_modes", "Poses per ligand", "9"),
    ("energy_range", "Energy range (kcal/mol)", "3"),
    ("seed", "Random seed", ""),
    ("granularity", "Grid spacing (Å)", "0.375"),
]

HELP = {
    "workflow": (
        "How to use DockingFlow",
        "Work through the tabs from left to right:\n\n"
        "1. Ligands: choose which ZINC22 molecules to screen, then click Create tranche list.\n"
        "2. Targets: choose your receptor and where on it to dock (the grid box).\n"
        "3. Docking: check the docking program and settings (the defaults are fine to start).\n"
        "4. Resources: how much of the server to use (pre-set to a recommendation).\n\n"
        "Then click Validate to check everything, and Run pipeline to start. You can close this window "
        "at any time; the run keeps going on the server. Reopen it with 'bash gui.sh'.\n\n"
        "Tip: start with the default 'Quick test' library to make sure everything works before a big screen."),
    "library": (
        "Ligand library",
        "ZINC22 is a free database of billions of purchasable compounds, split into 'tranches' by size "
        "(heavy atoms) and greasiness (logP), and further by charge.\n\n"
        "Pick a preset for a common kind of screen, or choose Custom and set the ranges yourself:\n\n"
        + "\n\n".join(f"- {name}: {p['help']}" for name, p in LIGAND_PRESETS.items())
        + "\n\nThe count under the ranges shows how many molecules you've selected. Docking takes very "
        "roughly 10-60 CPU-seconds per molecule, so millions of molecules need days to weeks even on a "
        "big server. Start small."),
    "hac": (
        "Heavy atoms",
        "The number of non-hydrogen atoms in a molecule: ZINC22's measure of size. Multiply by ~14 for a "
        "rough molecular weight (e.g. 25 heavy atoms is about 350 Da).\n\n"
        f"ZINC22 has 3D (dockable) structures for {HAC_MIN} to {HAC_MAX} heavy atoms. Typical ranges: "
        "fragments 8-19, lead-like 17-25, drug-like 17-29."),
    "logp": (
        "logP",
        "How greasy (lipophilic) a molecule is: the log of how it splits between oil and water. Negative = "
        "water-loving, higher = greasier.\n\n"
        "Most oral drugs have logP between about -1 and 5 (Lipinski's rule: <= 5). Lead-like libraries "
        "stay <= 3.5 to leave room for optimization."),
    "charges": (
        "Charges",
        "The molecule's net charge at physiological pH. ZINC22 files each protonation state separately.\n\n"
        "- Neutral only: simplest; good for a first screen.\n"
        "- Neutral and +/-1: covers most drugs (many are weak acids or bases).\n"
        "- All charges: includes highly charged molecules, which rarely make good drugs."),
    "tranche_list": (
        "Tranche list",
        "The list of ZINC22 files to download, made by 'Create tranche list' (or imported from a "
        "CartBlanche22 download file). You don't need to edit it."),
    "workdir": (
        "Work directory",
        "Where downloads, docking results and logs are stored. It can grow large (roughly 1-5 GB per "
        "million molecules), so choose a disk with plenty of space. Re-running with the same work "
        "directory continues where the last run stopped."),
    "targets": (
        "Receptor and grid box",
        "The receptor is your prepared protein structure (.pdbqt, e.g. made with ADFR's "
        "prepare_receptor or Meeko). The grid box is the region of it where ligands are docked.\n\n"
        "Known binding site: center the box on it, about 20-30 Å per side. The easiest way: 'Center on "
        "known site' with a ligand file from the same crystal structure (e.g. a co-crystallized "
        "inhibitor saved from PyMOL or ChimeraX).\n\n"
        "Unknown site (blind docking): 'Fit to whole protein' measures the receptor and covers all of it "
        "(each side capped at 80 Å). Slower and less precise, so raise exhaustiveness to 16-32 in the "
        "Docking tab.\n\n"
        "Center and size are in Ångströms, in the receptor file's own coordinates."),
    "blind": (
        "Fit to whole protein (blind docking)",
        "Reads every atom in the receptor file, centers the box on the protein and makes it big enough "
        "to cover it plus a 5 Å margin, up to 80 Å per side. Use this when you don't know where "
        "ligands bind."),
    "site": (
        "Center on known site",
        "Pick a file with a ligand already sitting in the binding site (PDB, PDBQT, MOL2 or SDF), e.g. the "
        "inhibitor from the crystal structure your receptor came from. The box is centered on it with a "
        "5 Å margin, at least 20 Å per side.\n\n"
        "The file must be in the same coordinate frame as the receptor, i.e. from the same structure."),
    "binary": (
        "Docking program",
        "VinaLC, the parallel version of AutoDock Vina used for docking. 'vinalc' works if it's on the "
        "server's PATH; otherwise Browse to it. The MPI launcher (usually 'mpirun') starts VinaLC on "
        "many cores; leave 'Skip MPI launcher' off unless you're testing without VinaLC."),
    "filter_percent": (
        "Keep top %",
        "After docking, each tranche's ligands are ranked by predicted binding strength and this top "
        "percentage is kept as hits (always at least one). All scores are still saved in ranked_all.tsv. "
        "Default 10."),
    "exhaustiveness": (
        "Exhaustiveness",
        "How hard Vina searches for each ligand's best pose. Higher is more reliable but slower (roughly "
        "proportional). 8 is the standard default; use 16-32 for large (blind-docking) boxes."),
    "num_modes": (
        "Poses per ligand",
        "How many alternative binding poses Vina reports per ligand. Only the best one is used for "
        "ranking. Default 9."),
    "energy_range": (
        "Energy range",
        "Only poses within this many kcal/mol of the best pose are reported. Doesn't affect ranking. "
        "Default 3."),
    "seed": (
        "Random seed",
        "Leave blank for a random search each run. Set a whole number to make runs reproducible."),
    "granularity": (
        "Grid spacing",
        "Spacing of Vina's precomputed energy grids, in Å. The default 0.375 is right for almost "
        "everyone; smaller is slower and uses much more memory."),
    "resources": (
        "CPU and memory",
        "How much of the server this run may use. The sliders start at a recommendation: physical cores "
        "minus what's already busy and a small reserve, and 80% of free memory. Lower them if you share "
        "the server.\n\n"
        "VinaLC runs one coordinator plus workers, each using about 'exhaustiveness' threads, so the "
        "number of workers is cores / exhaustiveness. Each worker needs memory for its grid box (bigger "
        "boxes need more); the panel below shows the resulting plan."),
    "settings_file": (
        "Settings file",
        "Where the GUI saves your docking settings, targets and CPU/memory choice (setup.txt), each time "
        "you click Validate or Run. The command-line pipeline reads the same file. You don't need to "
        "edit it by hand."),
}


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



class ScrollFrame(ttk.Frame):
    """A vertically scrollable area: put widgets in `.inner`. Keeps tabs usable in short windows."""

    def __init__(self, parent: tk.Widget, padding: int = 12) -> None:
        super().__init__(parent)
        self.canvas = tk.Canvas(self, highlightthickness=0, bd=0, background=P["card"])
        self.bar = ttk.Scrollbar(self, orient=tk.VERTICAL, command=self.canvas.yview)
        self.inner = ttk.Frame(self.canvas, padding=padding)
        self._win = self.canvas.create_window(0, 0, window=self.inner, anchor="nw")
        self.canvas.configure(yscrollcommand=self.bar.set)
        self.canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.bar.pack(side=tk.RIGHT, fill=tk.Y)
        self.inner.bind("<Configure>", lambda _e: self.canvas.configure(scrollregion=self.canvas.bbox("all")))
        self.canvas.bind("<Configure>", lambda e: self.canvas.itemconfigure(self._win, width=e.width))
        # Mouse wheel scrolls this area while the pointer is over it (X11 sends Button-4/5).
        self.bind("<Enter>", lambda _e: self._wheel(True))
        self.bind("<Leave>", lambda _e: self._wheel(False))

    def _wheel(self, on: bool) -> None:
        for seq in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
            if on:
                self.bind_all(seq, self._scroll)
            else:
                self.unbind_all(seq)

    def _scroll(self, e: tk.Event) -> None:
        if self.inner.winfo_reqheight() <= self.canvas.winfo_height():
            return
        step = -1 if (getattr(e, "num", 0) == 4 or getattr(e, "delta", 0) > 0) else 1
        self.canvas.yview_scroll(step, "units")


class DockingFlowApp(tk.Tk):
    def __init__(self, api: ApiClient) -> None:
        super().__init__()
        self.api = api
        self.title("DockingFlow")
        # Fit the screen (MobaXterm windows are often smaller than the X display suggests).
        w = min(1180, self.winfo_screenwidth() - 60)
        h = min(760, self.winfo_screenheight() - 100)
        self.geometry(f"{w}x{h}")
        self.minsize(760, 500)
        self._style()

        self.v = {name: tk.StringVar(self) for name in ("setup_path", "map_path", "workdir", "vinalc_bin", "mpirun_bin")}
        self.no_mpirun = tk.BooleanVar(self, value=False)

        # Ligands tab: preset + ranges + charges
        self.preset = tk.StringVar(self, value=DEFAULT_PRESET)
        self.hac_min, self.hac_max = tk.StringVar(self), tk.StringVar(self)
        self.logp_min, self.logp_max = tk.StringVar(self), tk.StringVar(self)
        self.charge_preset = tk.StringVar(self)
        self.charges = {c: tk.BooleanVar(self) for c in CHARGES}
        self._applying_preset = False
        self._zinc_job = None

        # Docking tab: settings saved into setup.txt
        self.settings = {k: tk.StringVar(self, value=d) for k, _, d in SETTINGS}
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
        self.apply_preset()
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
        # Pixel sizes (negative), not points: MobaXterm's X server often reports a high DPI, which
        # makes point-sized fonts render huge.
        base = tkfont.nametofont("TkDefaultFont")
        base.configure(size=-12)
        for name in ("TkTextFont", "TkMenuFont", "TkHeadingFont", "TkCaptionFont", "TkTooltipFont"):
            tkfont.nametofont(name).configure(family=base.cget("family"), size=-12)
        self.font_bold = base.copy()
        self.font_bold.configure(weight="bold")
        self.font_small = base.copy()
        self.font_small.configure(size=-11)
        self.font_title = base.copy()
        self.font_title.configure(size=-17, weight="bold")
        self.font_mono = tkfont.nametofont("TkFixedFont").copy()
        self.font_mono.configure(size=-11)
        self.option_add("*TCombobox*Listbox.font", base)  # the dropdown lists

        self.configure(background=P["bg"])
        st = ttk.Style(self)
        st.theme_use("clam")
        # Content widgets sit on white cards by default; the window frame is gray ("App.*").
        st.configure(".", background=P["card"], foreground=P["text"], bordercolor=P["border"],
                     lightcolor=P["card"], darkcolor=P["card"], troughcolor=P["bg"], focuscolor=P["accent"])
        st.configure("App.TFrame", background=P["bg"])
        st.configure("Muted.TLabel", foreground=P["muted"], font=self.font_small)
        st.configure("AppMuted.TLabel", background=P["bg"], foreground=P["muted"])
        st.configure("Title.TLabel", background=P["bg"], font=self.font_title)
        st.configure("Section.TLabel", font=self.font_bold)
        st.configure("Big.TLabel", font=self.font_bold, foreground=P["accent"])

        st.configure("TLabelframe", background=P["card"], bordercolor=P["border"], relief="solid", borderwidth=1)
        st.configure("TLabelframe.Label", background=P["card"], foreground=P["muted"], font=self.font_bold)
        st.configure("Card.TFrame", background=P["card"], bordercolor=P["border"], relief="solid", borderwidth=1)

        st.configure("TNotebook", background=P["bg"], borderwidth=0, tabmargins=(0, 0, 0, 0))
        st.configure("TNotebook.Tab", background=P["button"], foreground=P["muted"], padding=(12, 4),
                     bordercolor=P["border"])
        st.map("TNotebook.Tab", background=[("selected", P["card"])], foreground=[("selected", P["text"])],
               expand=[("selected", (0, 0, 0, 0))])

        st.configure("TButton", background=P["button"], padding=(8, 3), bordercolor=P["border"], relief="flat")
        st.map("TButton", background=[("active", P["button_hover"]), ("disabled", P["button"])],
               foreground=[("disabled", P["muted"])])
        st.configure("Accent.TButton", background=P["accent"], foreground="#ffffff", bordercolor=P["accent"],
                     font=self.font_bold)
        st.map("Accent.TButton", background=[("active", P["accent_hover"]), ("disabled", P["accent_disabled"])],
               foreground=[("disabled", "#ffffff")])
        st.configure("Danger.TButton", foreground=P["err"])
        st.map("Danger.TButton", background=[("active", P["err_bg"])])
        st.configure("Help.TButton", foreground=P["accent"], font=self.font_bold, padding=(4, 0), width=2,
                     background=P["info_bg"], bordercolor=P["info_bg"])
        st.map("Help.TButton", background=[("active", P["select"])])

        for widget in ("TEntry", "TSpinbox", "TCombobox"):
            st.configure(widget, fieldbackground="#ffffff", bordercolor=P["border"], padding=2,
                         lightcolor=P["border"], darkcolor=P["border"], arrowcolor=P["muted"])
            st.map(widget, bordercolor=[("focus", P["accent"])], lightcolor=[("focus", P["accent"])],
                   fieldbackground=[("readonly", "#ffffff")], background=[("readonly", P["button"])])
        st.configure("TCheckbutton", background=P["card"])
        st.map("TCheckbutton", background=[("active", P["card"])])

        st.configure("Treeview", background=P["card"], fieldbackground=P["card"], rowheight=20,
                     bordercolor=P["border"], borderwidth=0)
        st.map("Treeview", background=[("selected", P["select"])], foreground=[("selected", P["text"])])
        st.configure("Treeview.Heading", background=P["stripe"], foreground=P["muted"], font=self.font_bold,
                     relief="flat", padding=(4, 2))
        st.map("Treeview.Heading", background=[("active", P["button"])])
        st.configure("Accent.Horizontal.TProgressbar", background=P["accent"], troughcolor=P["neutral_bg"],
                     bordercolor=P["neutral_bg"], lightcolor=P["accent"], darkcolor=P["accent"], thickness=6)
        st.configure("TPanedwindow", background=P["bg"])
        for bar in ("Vertical.TScrollbar", "Horizontal.TScrollbar"):
            st.configure(bar, background=P["button"], troughcolor=P["card"], bordercolor=P["card"],
                         arrowcolor=P["muted"])

    # ---------- small building blocks ----------

    def help_button(self, parent: tk.Widget, key: str) -> ttk.Button:
        return ttk.Button(parent, text="?", style="Help.TButton", command=lambda: self.show_help(key))

    def show_help(self, key: str) -> None:
        title, text = HELP[key]
        win = tk.Toplevel(self)
        win.title(title)
        win.configure(background=P["card"])
        win.transient(self)
        win.resizable(False, False)
        frame = ttk.Frame(win, padding=14)
        frame.pack(fill=tk.BOTH, expand=True)
        ttk.Label(frame, text=title, style="Section.TLabel").pack(anchor="w", pady=(0, 6))
        ttk.Label(frame, text=text, wraplength=420, justify=tk.LEFT).pack(anchor="w")
        ok = ttk.Button(frame, text="Got it", style="Accent.TButton", command=win.destroy)
        ok.pack(anchor="e", pady=(12, 0))
        win.bind("<Escape>", lambda _e: win.destroy())
        win.bind("<Return>", lambda _e: win.destroy())
        ok.focus_set()

    def heading(self, parent: tk.Widget, text: str, help_key: Optional[str], pady=(10, 4)) -> ttk.Frame:
        row = ttk.Frame(parent)
        row.pack(fill=tk.X, pady=pady)
        ttk.Label(row, text=text, style="Section.TLabel").pack(side=tk.LEFT)
        if help_key:
            self.help_button(row, help_key).pack(side=tk.LEFT, padx=6)
        return row

    def path_row(self, parent: tk.Widget, var: tk.StringVar, kind: Optional[str]) -> ttk.Frame:
        row = ttk.Frame(parent)
        row.pack(fill=tk.X)
        ttk.Entry(row, textvariable=var).pack(side=tk.LEFT, fill=tk.X, expand=True)
        if kind:
            ttk.Button(row, text="Browse", command=lambda: self.browse(var, kind)).pack(side=tk.LEFT, padx=(4, 0))
        return row

    # ---------- layout ----------

    def _build(self) -> None:
        outer = ttk.Frame(self, style="App.TFrame", padding=(10, 8, 10, 10))
        outer.pack(fill=tk.BOTH, expand=True)

        header = ttk.Frame(outer, style="App.TFrame")
        header.pack(fill=tk.X, pady=(0, 8))
        ttk.Label(header, text="DockingFlow", style="Title.TLabel").pack(side=tk.LEFT)
        self.subtitle = ttk.Label(header, text="virtual screening with VinaLC", style="AppMuted.TLabel")
        self.subtitle.pack(side=tk.LEFT, padx=(8, 0), pady=(4, 0))
        self.pill = tk.Label(header, text="Idle", font=self.font_bold, padx=10, pady=2)
        self.pill.pack(side=tk.RIGHT)
        ttk.Button(header, text="How to use", command=lambda: self.show_help("workflow")).pack(side=tk.RIGHT, padx=8)
        self.set_phase("idle")

        panes = ttk.PanedWindow(outer, orient=tk.HORIZONTAL)
        panes.pack(fill=tk.BOTH, expand=True)
        left = ttk.Frame(panes, style="App.TFrame", width=400)
        right = ttk.Frame(panes, style="App.TFrame")
        panes.add(left, weight=0)
        panes.add(right, weight=1)

        tabs = ttk.Notebook(left)
        tabs.pack(fill=tk.BOTH, expand=True, padx=(0, 10))
        for title, build in (("1 Ligands", self._build_ligands), ("2 Targets", self._build_targets_tab),
                             ("3 Docking", self._build_docking), ("4 Resources", self._build_resources)):
            area = ScrollFrame(tabs)
            build(area.inner)
            tabs.add(area, text=title)

        actions = ttk.Frame(left, style="Card.TFrame", padding=8)
        actions.pack(fill=tk.X, padx=(0, 10), pady=(8, 0))
        row = ttk.Frame(actions)
        row.pack(fill=tk.X)
        ttk.Button(row, text="Validate", command=self.validate).pack(side=tk.LEFT)
        self.run_btn = ttk.Button(row, text="Run pipeline", style="Accent.TButton", command=self.start_run)
        self.run_btn.pack(side=tk.LEFT, padx=6)
        ttk.Button(row, text="Nuke", style="Danger.TButton", command=self.nuke).pack(side=tk.RIGHT)
        ttk.Button(row, text="Clean", command=self.clean).pack(side=tk.RIGHT, padx=6)
        self.message = tk.Label(actions, text="Loading...", wraplength=360, justify=tk.LEFT, anchor="w",
                                padx=8, pady=6, font=self.font_small)
        self.message.pack(fill=tk.X, pady=(8, 0))
        self.set_message("Loading...")

        self._build_status(right)

    def set_phase(self, phase: str) -> None:
        text, fg, bg = PHASE_PILL.get(phase, PHASE_PILL["idle"])
        self.pill.configure(text=text, fg=fg, bg=bg)

    def _build_ligands(self, f: ttk.Frame) -> None:
        self.heading(f, "Ligand library", "library", pady=(0, 4))
        combo = ttk.Combobox(f, textvariable=self.preset, state="readonly",
                             values=list(LIGAND_PRESETS) + [CUSTOM])
        combo.pack(fill=tk.X)
        combo.bind("<<ComboboxSelected>>", lambda _e: self.apply_preset())
        self.preset_note = ttk.Label(f, text="", style="Muted.TLabel", wraplength=340, justify=tk.LEFT)
        self.preset_note.pack(fill=tk.X, pady=(3, 0))

        grid = ttk.Frame(f)
        grid.pack(fill=tk.X, pady=(8, 0))
        hac_values = [str(n) for n in range(HAC_MIN, HAC_MAX + 1)]
        logp_values = [f"{v:g}" for v in LOGP_BINS]
        self.mw_note = ttk.Label(grid, text="", style="Muted.TLabel")
        for r, (label, lo, hi, values, key) in enumerate((
            ("Heavy atoms", self.hac_min, self.hac_max, hac_values, "hac"),
            ("logP", self.logp_min, self.logp_max, logp_values, "logp"),
        )):
            ttk.Label(grid, text=label).grid(row=r * 2, column=0, sticky="w", pady=(2, 0))
            for col, var in ((1, lo), (3, hi)):
                cb = ttk.Combobox(grid, textvariable=var, values=values, state="readonly", width=6)
                cb.grid(row=r * 2, column=col, padx=2, pady=(2, 0))
                cb.bind("<<ComboboxSelected>>", lambda _e: self.ranges_changed())
            ttk.Label(grid, text="to").grid(row=r * 2, column=2)
            self.help_button(grid, key).grid(row=r * 2, column=4, padx=(6, 0))
        self.mw_note.grid(row=1, column=1, columnspan=4, sticky="w")
        ttk.Label(grid, text="oral drugs: about -1 to 5", style="Muted.TLabel").grid(
            row=3, column=1, columnspan=4, sticky="w")

        ttk.Label(grid, text="Charges").grid(row=4, column=0, sticky="w", pady=(6, 0))
        ccombo = ttk.Combobox(grid, textvariable=self.charge_preset, state="readonly",
                              values=list(CHARGE_PRESETS) + [CUSTOM], width=20)
        ccombo.grid(row=4, column=1, columnspan=3, sticky="w", padx=2, pady=(6, 0))
        ccombo.bind("<<ComboboxSelected>>", lambda _e: self.apply_charge_preset())
        self.help_button(grid, "charges").grid(row=4, column=4, padx=(6, 0), pady=(6, 0))
        boxes = ttk.Frame(f)
        boxes.pack(fill=tk.X, pady=(4, 0))
        for c, label in CHARGES.items():
            ttk.Checkbutton(boxes, text=label, variable=self.charges[c], command=self.charges_changed).pack(side=tk.LEFT)

        summary = ttk.Frame(f, style="Card.TFrame", padding=8)
        summary.pack(fill=tk.X, pady=(10, 0))
        self.zinc_count = ttk.Label(summary, text="", style="Big.TLabel")
        self.zinc_count.pack(anchor="w")
        self.zinc_label = ttk.Label(summary, text="", style="Muted.TLabel", wraplength=330, justify=tk.LEFT)
        self.zinc_label.pack(anchor="w")
        ttk.Button(summary, text="Create tranche list", style="Accent.TButton", command=self.zinc_build).pack(
            anchor="w", pady=(6, 0))

        self.heading(f, "Tranche list", "tranche_list")
        self.path_row(f, self.v["map_path"], "file")
        ttk.Button(f, text="Import a CartBlanche22 file instead...", command=self.import_zinc).pack(anchor="w", pady=(4, 0))
        self.heading(f, "Work directory", "workdir")
        self.path_row(f, self.v["workdir"], "folder")

    def _build_targets_tab(self, f: ttk.Frame) -> None:
        self.heading(f, "Receptors and grid boxes", "targets", pady=(0, 4))
        ttk.Label(f, style="Muted.TLabel", wraplength=340, justify=tk.LEFT, text=(
            "Known binding site: box of ~20-30 Å around it ('Center on known site').\n"
            "Unknown site: blind docking over the whole protein, at most 80 x 80 x 80 Å "
            "('Fit to whole protein').")).pack(fill=tk.X)
        self.targets_frame = ttk.Frame(f)
        self.targets_frame.pack(fill=tk.X, pady=6)
        buttons = ttk.Frame(f)
        buttons.pack(fill=tk.X)
        ttk.Button(buttons, text="+ Add receptor", command=self.add_target).pack(side=tk.LEFT)
        ttk.Button(buttons, text="Save targets", command=self.save_targets).pack(side=tk.LEFT, padx=6)

    def _build_docking(self, f: ttk.Frame) -> None:
        self.heading(f, "Docking program", "binary", pady=(0, 4))
        ttk.Label(f, text="VinaLC (name on PATH, or full path)", style="Muted.TLabel").pack(anchor="w")
        self.path_row(f, self.v["vinalc_bin"], "file")
        ttk.Label(f, text="MPI launcher", style="Muted.TLabel").pack(anchor="w", pady=(4, 0))
        self.path_row(f, self.v["mpirun_bin"], None)
        ttk.Checkbutton(f, text="Skip MPI launcher (testing only)", variable=self.no_mpirun).pack(anchor="w", pady=(4, 0))

        self.heading(f, "Docking settings", None)
        ttk.Label(f, text="Defaults are fine for most screens. Saved when you Validate or Run.",
                  style="Muted.TLabel", wraplength=340).pack(anchor="w")
        grid = ttk.Frame(f)
        grid.pack(fill=tk.X, pady=(4, 0))
        for r, (key, label, default) in enumerate(SETTINGS):
            ttk.Label(grid, text=label).grid(row=r, column=0, sticky="w", pady=2)
            ttk.Entry(grid, textvariable=self.settings[key], width=8, justify=tk.RIGHT).grid(row=r, column=1, padx=6)
            self.help_button(grid, key).grid(row=r, column=2)
            if default:
                ttk.Label(grid, text=f"default {default}", style="Muted.TLabel").grid(row=r, column=3, sticky="w", padx=6)
        ttk.Button(f, text="Reset to defaults", command=self.reset_settings).pack(anchor="w", pady=(6, 0))
        self.settings["exhaustiveness"].trace_add("write", lambda *_: self.schedule_plan())

        self.heading(f, "Settings file", "settings_file")
        self.path_row(f, self.v["setup_path"], "file")
        self.v["setup_path"].trace_add("write", lambda *_: self.after_idle(self._setup_changed))

    def _build_resources(self, f: ttk.Frame) -> None:
        self.heading(f, "CPU and memory", "resources", pady=(0, 4))
        buttons = ttk.Frame(f)
        buttons.pack(fill=tk.X)
        ttk.Button(buttons, text="Detect this machine", command=self.detect_resources).pack(side=tk.LEFT)
        ttk.Button(buttons, text="Load server report...", command=self.load_report).pack(side=tk.LEFT, padx=6)

        self.specs = ttk.Label(f, text="Detecting...", justify=tk.LEFT, padding=(0, 6), style="Muted.TLabel")
        self.specs.pack(fill=tk.X)

        ttk.Label(f, text="CPU cores to use").pack(anchor="w")
        self.cores_scale = tk.Scale(f, from_=1, to=1, orient=tk.HORIZONTAL, variable=self.cores,
                                    command=lambda _: self.schedule_plan(), **self._scale_look())
        self.cores_scale.pack(fill=tk.X)
        ttk.Label(f, text="Memory budget (GB)").pack(anchor="w", pady=(4, 0))
        self.mem_scale = tk.Scale(f, from_=1, to=1, orient=tk.HORIZONTAL, variable=self.memory_gb,
                                  command=lambda _: self.schedule_plan(), **self._scale_look())
        self.mem_scale.pack(fill=tk.X)

        plan = ttk.Frame(f, style="Card.TFrame", padding=8)
        plan.pack(fill=tk.X, pady=(8, 0))
        self.plan_label = ttk.Label(plan, text="", justify=tk.LEFT, wraplength=330)
        self.plan_label.pack(fill=tk.X)
        ttk.Button(f, text="Use recommended", command=self.use_recommended).pack(anchor="w", pady=(6, 0))

    def _scale_look(self) -> Dict[str, Any]:
        return dict(highlightthickness=0, bg=P["card"], fg=P["text"], troughcolor=P["neutral_bg"],
                    activebackground=P["accent_hover"], sliderrelief=tk.FLAT, bd=0, sliderlength=18, width=10,
                    font=self.font_small)

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
        card = ttk.Frame(parent, style="Card.TFrame", padding=(10, 8))
        card.pack(**pack)
        head = ttk.Frame(card)
        head.pack(fill=tk.X, pady=(0, 4))
        ttk.Label(head, text=title, style="Section.TLabel").pack(side=tk.LEFT)
        card.head = head
        return card

    def _build_status(self, right: ttk.Frame) -> None:
        box = self._card(right, "Tranches", fill=tk.X)
        ttk.Button(box.head, text="Details...", command=self.show_details).pack(side=tk.RIGHT)
        self.progress_text = ttk.Label(box.head, text="", style="Muted.TLabel")
        self.progress_text.pack(side=tk.RIGHT, padx=8)
        self.progress = ttk.Progressbar(box, style="Accent.Horizontal.TProgressbar", mode="determinate")
        self.progress.pack(fill=tk.X, pady=(0, 6))
        self.tranche_tree = self._tree(
            box, [("Tranche", 90), ("LogP", 50, "e"), ("Size", 70), ("Status", 110), ("Error", 360)], 5)
        self.tranche_tree.tag_configure("DONE", foreground=P["ok"])
        self.tranche_tree.tag_configure("FAILED", foreground=P["err"])
        self.tranche_tree.tag_configure("ACTIVE", foreground=P["info"])
        self.tranche_tree.bind("<Double-1>", lambda _e: self.show_details())
        ttk.Label(box, text="Select a tranche and click Details (or double-click it) to see its logs and errors.",
                  style="Muted.TLabel").pack(anchor="w", pady=(4, 0))

        box = self._card(right, "Log", fill=tk.BOTH, pady=8)
        frame = ttk.Frame(box)
        frame.pack(fill=tk.BOTH, expand=True)
        self.log = tk.Text(frame, height=7, wrap=tk.WORD, state=tk.DISABLED, font=self.font_mono,
                           bg=P["log_bg"], fg=P["log_fg"], relief=tk.FLAT, padx=8, pady=6,
                           highlightthickness=0, insertbackground=P["log_fg"])
        scroll = ttk.Scrollbar(frame, orient=tk.VERTICAL, command=self.log.yview)
        self.log.configure(yscrollcommand=scroll.set)
        self.log.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)

        box = self._card(right, "Top hits", fill=tk.BOTH, expand=True)
        ttk.Label(box.head, text="strongest predicted binders first (more negative = stronger)",
                  style="Muted.TLabel").pack(side=tk.LEFT, padx=8)
        self.hits_tree = self._tree(
            box, [("Ligand", 220), ("Affinity (kcal/mol)", 120, "e"), ("Receptor", 130), ("Tranche", 90)], 8)


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
                self.set_message("New here? Click 'How to use' (top right), then work through tabs 1-4 "
                                 "and click Validate. Every '?' explains that setting.")
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
        for key, _, default in SETTINGS:
            self.settings[key].set(saved.get(key) or default)
        self._saved_budget = (saved.get("cores"), saved.get("memory_gb"))
        self.apply_saved_budget()

    def reset_settings(self) -> None:
        for key, _, default in SETTINGS:
            self.settings[key].set(default)
        self.set_message("Docking settings reset to defaults (saved on Validate/Run).", "info")

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
            self.set_message(f"Docking settings (tab 3): {res['message']}", "err")
        return res["ok"]

    def save_all(self) -> bool:
        """Save targets (if edited) and settings; False if anything was invalid."""
        if self.targets_dirty and not self.save_targets():
            return False
        return self.save_settings()

    # ---------- ligand library: presets + ranges + charges ----------

    def apply_preset(self) -> None:
        name = self.preset.get()
        preset = LIGAND_PRESETS.get(name)
        if not preset:  # Custom: keep the current ranges
            self.preset_note.configure(text="Set the ranges and charges below yourself.")
            return
        self._applying_preset = True
        self.hac_min.set(str(preset["hac"][0]))
        self.hac_max.set(str(preset["hac"][1]))
        self.logp_min.set(f"{preset['logp'][0]:g}")
        self.logp_max.set(f"{preset['logp'][1]:g}")
        self.charge_preset.set(preset["charges"])
        self.apply_charge_preset()
        self._applying_preset = False
        self.preset_note.configure(text=preset["help"])
        self.ranges_changed(from_preset=True)

    def apply_charge_preset(self) -> None:
        wanted = CHARGE_PRESETS.get(self.charge_preset.get())
        if wanted is None:
            return
        for c, var in self.charges.items():
            var.set(c in wanted)
        self.charges_changed(from_preset=True)

    def charges_changed(self, from_preset: bool = False) -> None:
        chosen = sorted(c for c, v in self.charges.items() if v.get())
        match = [name for name, cs in CHARGE_PRESETS.items() if sorted(cs) == chosen]
        self.charge_preset.set(match[0] if match else CUSTOM)
        if not from_preset:
            self.ranges_changed()

    def ranges_changed(self, from_preset: bool = False) -> None:
        if not from_preset and not self._applying_preset:
            self.preset.set(CUSTOM)
            self.preset_note.configure(text="Custom selection.")
        try:
            lo, hi = int(self.hac_min.get()), int(self.hac_max.get())
            mw = f"{lo * DA_PER_HEAVY_ATOM}" + (f"-{hi * DA_PER_HEAVY_ATOM}" if hi != lo else "")
            self.mw_note.configure(text=f"about {mw} Da")
        except ValueError:
            pass
        self.schedule_zinc_summary()

    def _zinc_filters(self) -> Optional[Dict[str, Any]]:
        try:
            f = {
                "hac_min": int(self.hac_min.get()), "hac_max": int(self.hac_max.get()),
                "logp_min": float(self.logp_min.get()), "logp_max": float(self.logp_max.get()),
                "charges": [c for c, v in self.charges.items() if v.get()],
            }
        except ValueError:
            return None
        if f["hac_min"] > f["hac_max"] or f["logp_min"] > f["logp_max"] or not f["charges"]:
            return None
        return f

    def schedule_zinc_summary(self) -> None:
        if self._zinc_job:
            self.after_cancel(self._zinc_job)
        self._zinc_job = self.after(300, self.refresh_zinc_summary)

    def refresh_zinc_summary(self) -> None:
        self._zinc_job = None
        filters = self._zinc_filters()
        if filters is None:
            self.zinc_count.configure(text="Nothing selected")
            self.zinc_label.configure(text="Each 'from' must be <= its 'to', and pick at least one charge.")
            return
        self.zinc_count.configure(text="Counting...")
        self.zinc_label.configure(text="(the first time downloads ZINC22's index, ~10 MB)")
        self.update_idletasks()
        res = self.call("zinc_summary", filters)
        if res is None:
            return
        if not res["ok"]:
            self.zinc_count.configure(text="Couldn't count")
            self.zinc_label.configure(text=res["message"])
            return
        n = res["molecules"]
        self.zinc_count.configure(text=f"~{n:,} molecules")
        note = f"in {res['tranches']:,} tranche(s)."
        if n >= LARGE_SCREEN:
            note += " That's a very large screen (days to weeks); consider narrowing it or starting with a quick test."
        self.zinc_label.configure(text=note)

    def zinc_build(self) -> None:
        filters = self._zinc_filters()
        if filters is None:
            self.set_message("Check the ligand ranges: each 'from' must be <= its 'to', with at least one charge.", "err")
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
        self.set_message(res["message"] + (" Next: tab 2, Targets." if res["ok"] else ""), "ok" if res["ok"] else "err")

    # ---------- file pickers + ZINC import ----------

    def browse(self, var: tk.StringVar, kind: str, title: Optional[str] = None,
               filetypes: Optional[list] = None) -> Optional[str]:
        start = Path(var.get() or self.v["setup_path"].get() or ".").expanduser()
        initial = str(start if start.is_dir() else start.parent)
        if kind == "folder":
            path = filedialog.askdirectory(initialdir=initial, title=title or "Choose a folder")
        else:
            path = filedialog.askopenfilename(initialdir=initial, title=title or "Choose a file",
                                              filetypes=filetypes or [("All files", "*")])
        if path:
            var.set(path)
        return path or None

    def import_zinc(self) -> None:
        initial = str(Path(self.v["map_path"].get() or ".").expanduser().parent)
        path = filedialog.askopenfilename(
            initialdir=initial, title="ZINC22 download file from CartBlanche22",
            filetypes=[("Download scripts", "*.curl *.wget *.aws *.txt"), ("All files", "*")])
        if not path:
            return
        res = self.call("import_zinc_downloader", path)
        if res is None:
            return
        if res["ok"]:
            self.v["map_path"].set(res["map_path"])
            self.preset.set(CUSTOM)
            self.preset_note.configure(text="Using the imported file's tranches (the ranges above are ignored).")
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
        t = t or {"receptor": "", "center": ["", "", ""], "size": ["22", "22", "22"]}
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
            top = ttk.Frame(box)
            top.pack(fill=tk.X)
            ttk.Entry(top, textvariable=row["receptor"]).pack(side=tk.LEFT, fill=tk.X, expand=True)
            ttk.Button(top, text="Browse", command=lambda v=row["receptor"]: self.browse(
                v, "file", "Prepared receptor (.pdbqt)", [("PDBQT", "*.pdbqt"), ("All files", "*")])).pack(
                side=tk.LEFT, padx=(4, 0))
            ttk.Button(top, text="Remove", style="Danger.TButton", command=lambda k=i: self.remove_target(k)).pack(
                side=tk.LEFT, padx=(4, 0))

            grid = ttk.Frame(box)
            grid.pack(fill=tk.X, pady=(4, 0))
            for c in range(1, 4):
                grid.columnconfigure(c, weight=1)
            for k, axis in enumerate("xyz"):
                ttk.Label(grid, text=axis, style="Muted.TLabel", anchor="center").grid(row=0, column=1 + k, sticky="ew")
            for r, (label, key) in enumerate((("Center", "center"), ("Size (Å)", "size")), start=1):
                ttk.Label(grid, text=label).grid(row=r, column=0, sticky="w", padx=(0, 4))
                for k, var in enumerate(row[key]):
                    ttk.Entry(grid, textvariable=var, width=8, justify=tk.RIGHT).grid(
                        row=r, column=1 + k, sticky="ew", padx=1, pady=1)

            fit = ttk.Frame(box)
            fit.pack(fill=tk.X, pady=(6, 0))
            ttk.Button(fit, text="Fit to whole protein", command=lambda r=row: self.fit_blind(r)).pack(side=tk.LEFT)
            self.help_button(fit, "blind").pack(side=tk.LEFT, padx=(3, 8))
            ttk.Button(fit, text="Center on known site...", command=lambda r=row: self.fit_site(r)).pack(side=tk.LEFT)
            self.help_button(fit, "site").pack(side=tk.LEFT, padx=3)

    def _apply_box(self, row: Dict[str, Any], res: Optional[Dict[str, Any]]) -> None:
        if res is None:
            return
        if not res["ok"]:
            self.set_message(res["message"], "err")
            return
        for var, v in zip(row["center"], res["center"]):
            var.set(f"{v:g}")
        for var, v in zip(row["size"], res["size"]):
            var.set(f"{v:g}")
        self.set_message(res["message"] + " Saved when you Validate or Run.", "info" if res.get("capped") else "ok")
        self.schedule_plan()

    def fit_blind(self, row: Dict[str, Any]) -> None:
        receptor = row["receptor"].get().strip()
        if not receptor:
            self.set_message("Choose the receptor file first (Browse), then fit the box to it.", "err")
            return
        self._apply_box(row, self.call("fit_box_to_receptor", receptor))

    def fit_site(self, row: Dict[str, Any]) -> None:
        start = row["receptor"].get() or self.v["setup_path"].get() or "."
        path = filedialog.askopenfilename(
            initialdir=str(Path(start).expanduser().parent), title="Ligand in the binding site (same structure as the receptor)",
            filetypes=[("Structures", "*.pdb *.pdbqt *.mol2 *.sdf *.mol *.xyz"), ("All files", "*")])
        if path:
            self._apply_box(row, self.call("fit_box_to_site", path))

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
        if not res["ok"]:
            self.set_message(f"Targets (tab 2): {res['message']}", "err")
        else:
            self.set_message(res["message"], "info" if res["warnings"] else "ok")
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
            f"VinaLC will run mpirun -np {p['mpi_ranks']}: 1 coordinator + {p['workers']} workers "
            f"x ~{p['exhaustiveness']} threads.",
            f"Estimated memory: {p['est_memory_gb']:.1f} GB (~{p['rank_memory_gb']:.1f} GB per worker for your grid box).",
            f"Limited by: {p['limited_by']}.",
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
            self.set_message(res["message"] + ("\nLooks good: click Run pipeline." if res["ok"] else ""),
                             "ok" if res["ok"] else "err")

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
