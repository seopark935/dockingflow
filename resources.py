"""Machine resource detection and CPU/memory recommendations for a docking run.

Used by the GUI's Resources panel. Resources come from one of two places:

  - `detect_local()`: the machine this process is running on (the right
    choice when the GUI itself runs on the docking server).
  - `parse_server_report(text)`: the `DF_*` summary block printed by
    `server_check.sh`, so a run can be sized for a remote server from a
    laptop and the result written to setup.txt for the CLI.

`recommend()` turns those into a suggested core/memory budget, and
`plan()` shows what a chosen budget means for VinaLC (via
`io_parse.plan_mpi_ranks`).
"""
from __future__ import annotations

import os
import platform
import re
import subprocess
from pathlib import Path
from typing import Any

import io_parse


def _sysctl(name: str) -> str | None:
    try:
        return subprocess.run(["sysctl", "-n", name], capture_output=True, text=True, timeout=5).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return None


def _linux_physical_cores() -> int | None:
    """Count distinct (physical id, core id) pairs in /proc/cpuinfo."""
    try:
        text = Path("/proc/cpuinfo").read_text()
    except OSError:
        return None
    pairs = set()
    phys = core = None
    for line in text.splitlines() + [""]:
        if not line.strip():
            if core is not None:
                pairs.add((phys, core))
            phys = core = None
        elif line.startswith("physical id"):
            phys = line.split(":", 1)[1].strip()
        elif line.startswith("core id"):
            core = line.split(":", 1)[1].strip()
    return len(pairs) or None


def _linux_meminfo_kb() -> dict[str, int]:
    try:
        text = Path("/proc/meminfo").read_text()
    except OSError:
        return {}
    out = {}
    for line in text.splitlines():
        key, _, rest = line.partition(":")
        if rest.strip().endswith("kB"):
            out[key] = int(rest.split()[0])
    return out


def detect_local() -> dict[str, Any]:
    """Resources of the machine this process runs on, in the same shape as `parse_server_report`."""
    try:
        logical = len(os.sched_getaffinity(0))  # respects taskset/cgroup CPU limits
    except AttributeError:
        logical = os.cpu_count() or 1

    physical = total_gb = available_gb = None
    if platform.system() == "Linux":
        physical = _linux_physical_cores()
        mem = _linux_meminfo_kb()
        if "MemTotal" in mem:
            total_gb = mem["MemTotal"] / 1024**2
        if "MemAvailable" in mem:
            available_gb = mem["MemAvailable"] / 1024**2
    elif platform.system() == "Darwin":
        physical = int(_sysctl("hw.physicalcpu") or 0) or None
        memsize = _sysctl("hw.memsize")
        total_gb = int(memsize) / 1024**3 if memsize else None

    try:
        load1 = os.getloadavg()[0]
    except OSError:
        load1 = None

    return {
        "source": "this machine",
        "hostname": platform.node(),
        "logical_cpus": logical,
        "physical_cores": physical or logical,
        "mem_total_gb": total_gb,
        "mem_available_gb": available_gb,
        "load1": load1,
    }


def parse_server_report(text: str) -> dict[str, Any]:
    """Parse the `DF_KEY=VALUE` summary block from `server_check.sh` output.

    Raises:
        ValueError: If the text has no usable CPU count (e.g. it isn't a
            server_check.sh report, or the summary block was cut off).
    """
    vals = dict(re.findall(r"^DF_([A-Z0-9_]+)=(\S*)\s*$", text, flags=re.M))

    def num(key: str) -> float | None:
        try:
            return float(vals[key])
        except (KeyError, ValueError):
            return None

    logical = num("LOGICAL_CPUS")
    if not logical:
        raise ValueError("No DF_LOGICAL_CPUS line found. Is this the output of server_check.sh?")
    physical = num("PHYSICAL_CORES")
    total_kb, avail_kb = num("MEM_TOTAL_KB"), num("MEM_AVAILABLE_KB")
    disk_kb = num("DISK_FREE_KB")
    return {
        "source": "server report",
        "hostname": vals.get("HOSTNAME", "?"),
        "logical_cpus": int(logical),
        "physical_cores": int(physical) if physical else int(logical),
        "mem_total_gb": total_kb / 1024**2 if total_kb else None,
        "mem_available_gb": avail_kb / 1024**2 if avail_kb else None,
        "load1": num("LOAD1"),
        "disk_free_gb": disk_kb / 1024**2 if disk_kb else None,
        "inodes_free": int(num("INODES_FREE")) if num("INODES_FREE") else None,
    }


def recommend(res: dict[str, Any]) -> dict[str, Any]:
    """Suggest a core and memory budget that leaves the machine usable.

    Cores: physical cores (Vina's search threads gain little from
    hyperthreads), minus whatever is already busy (1-minute load average),
    minus a small reserve (~5%, at least 1) for the OS, MPI master, and I/O.
    Memory: 80% of currently available RAM (or of total, if availability
    is unknown).
    """
    physical = res["physical_cores"]
    busy = int(round(res.get("load1") or 0))
    reserve = max(1, physical // 20)
    cores = max(1, physical - busy - reserve)

    mem_basis = res.get("mem_available_gb") or res.get("mem_total_gb")
    memory_gb = round(mem_basis * 0.8, 1) if mem_basis else None

    notes = []
    if busy:
        notes.append(f"~{busy} core(s) already busy (load average {res['load1']:.1f}); left free.")
    if res["logical_cpus"] > physical:
        notes.append(f"{res['logical_cpus']} logical CPUs are hyperthreads of {physical} physical cores; sized for physical.")
    return {"cores": cores, "memory_gb": memory_gb, "notes": notes}


def plan(setup_path: str, cores: int, memory_gb: float | None) -> dict[str, Any]:
    """What a given core/memory budget means for VinaLC, using setup.txt's boxes and settings."""
    setup = io_parse.load_setup(Path(setup_path))
    targets = io_parse.load_docking_targets(setup)
    exhaustiveness = int(setup.get("exhaustiveness", io_parse.DEFAULT_EXHAUSTIVENESS))
    granularity = float(setup.get("granularity", io_parse.DEFAULT_GRANULARITY))

    rank_gb = io_parse.estimate_rank_memory_gb(targets, granularity)
    ranks = io_parse.plan_mpi_ranks(cores, exhaustiveness, memory_gb, rank_gb)
    workers = ranks - 1
    cpu_workers = max(1, cores // exhaustiveness)
    return {
        "mpi_ranks": ranks,
        "workers": workers,
        "exhaustiveness": exhaustiveness,
        "rank_memory_gb": rank_gb,
        "est_memory_gb": io_parse.RANK_OVERHEAD_GB + workers * rank_gb,
        "limited_by": "memory" if workers < cpu_workers else "cores",
    }


def set_setup_keys(setup_path: str, updates: dict[str, str | None]) -> None:
    """Set (or, for a `None` value, remove) `KEY=VALUE` lines in setup.txt in place.

    Every other line, including comments, is kept as-is; keys not already
    present are appended.
    """
    path = Path(setup_path).expanduser()
    pending = dict(updates)
    out = []
    for line in path.read_text().splitlines():
        key = line.split("=", 1)[0].strip()
        if "=" in line and not line.lstrip().startswith("#") and key in pending:
            val = pending.pop(key)
            if val is not None:
                out.append(f"{key}={val}")
        else:
            out.append(line)
    out += [f"{k}={v}" for k, v in pending.items() if v is not None]
    path.write_text("\n".join(out) + "\n")


def write_budget_to_setup(setup_path: str, cores: int, memory_gb: float | None) -> None:
    """Set `cores=` and `memory_gb=` in setup.txt (removing `memory_gb` when there's no budget)."""
    set_setup_keys(setup_path, {
        "cores": str(cores),
        "memory_gb": f"{memory_gb:g}" if memory_gb is not None else None,
    })
