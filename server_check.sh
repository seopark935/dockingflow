#!/bin/bash
# Report the hardware/software facts dockingflow depends on.
#
# Run on the docking server (from the directory you'll run the pipeline in):
#   bash server_check.sh 2>&1 | tee server_report.txt
# Optional: pass the planned workdir to check its disk space/inodes:
#   bash server_check.sh /path/to/workdir 2>&1 | tee server_report.txt
#
# Read-only apart from a 2-process `mpirun hostname` smoke test.

WORKDIR="${1:-$PWD}"
section() { printf '\n===== %s =====\n' "$1"; }
have() { command -v "$1" >/dev/null 2>&1; }

section "Host / OS"
hostname
uname -a
cat /etc/os-release 2>/dev/null | head -4
uptime

section "CPU (cores / threads / NUMA)"
echo "nproc (usable by this shell): $(nproc 2>/dev/null)"
if have lscpu; then
    lscpu | grep -E '^(Architecture|CPU\(s\)|On-line|Thread\(s\) per core|Core\(s\) per socket|Socket\(s\)|NUMA node\(s\)|Model name|CPU max MHz|L3 cache)'
else
    sysctl -n machdep.cpu.brand_string hw.physicalcpu hw.logicalcpu 2>/dev/null
fi

section "Memory"
free -h 2>/dev/null || vm_stat
echo "ulimit -n (open files): $(ulimit -n)   ulimit -u (processes/threads): $(ulimit -u)"

section "Current load / other users (is the machine shared?)"
who | awk '{print $1}' | sort | uniq -c
ps -eo user,pcpu,pmem,comm --sort=-pcpu 2>/dev/null | head -8

section "Disk space + inodes for workdir: $WORKDIR"
df -h "$WORKDIR" 2>/dev/null
df -i "$WORKDIR" 2>/dev/null
df -h /tmp 2>/dev/null | tail -1

section "Job scheduler (empty = none)"
for s in sbatch squeue qsub bsub; do have "$s" && echo "found: $s -> $(command -v $s)"; done
[ -n "$SLURM_JOB_ID$PBS_JOBID" ] && echo "inside a job: SLURM=$SLURM_JOB_ID PBS=$PBS_JOBID"

section "Environment modules"
if type module >/dev/null 2>&1; then
    module list 2>&1 | head -20
    module avail 2>&1 | grep -iE 'mpi|boost|python|vina' | head -20
else
    echo "no 'module' command"
fi

section "Python"
for py in python3 python; do have "$py" && echo "$py -> $(command -v $py): $($py --version 2>&1)"; done

section "Shell tools used by the pipeline"
for t in bash curl tar gzip; do
    have "$t" && echo "$t -> $(command -v $t): $($t --version 2>&1 | head -1)" || echo "MISSING: $t"
done

section "VinaLC"
VINALC="$(command -v vinalc || true)"
if [ -n "$VINALC" ]; then
    echo "vinalc -> $VINALC"
    have ldd && ldd "$VINALC" | grep -iE 'mpi|boost|not found'
else
    echo "vinalc not on PATH (note its full path and pass it via --vinalc-bin)"
fi

section "MPI"
if have mpirun; then
    echo "mpirun -> $(command -v mpirun)"
    mpirun --version 2>&1 | head -3
    echo "--- smoke test: mpirun -np 2 hostname"
    timeout 60 mpirun -np 2 hostname 2>&1 | head -5
else
    echo "mpirun not on PATH"
fi
[ "$(id -u)" = 0 ] && echo "NOTE: running as root -- Open MPI needs --allow-run-as-root"

section "GUI display (bash gui.sh)"
if [ -n "$DISPLAY" ]; then
    echo "X11 forwarding: ON (DISPLAY=$DISPLAY) -- gui.sh can open a window on your screen"
else
    echo "X11 forwarding: OFF -- in MobaXterm enable SSH > Advanced SSH settings > X11-Forwarding"
fi
found_browser=""
for b in chromium chromium-browser google-chrome google-chrome-stable firefox; do
    have "$b" && { echo "browser: $b -> $(command -v $b)"; found_browser=1; }
done
[ -z "$found_browser" ] && echo "browser: none (run 'bash gui.sh setup' for a small built-in window)"
[ -x .venv/bin/python3 ] && .venv/bin/python3 -c "import webview" 2>/dev/null \
    && echo "pywebview window: installed (.venv)"

section "Network access to ZINC"
curl -s -m 20 -o /dev/null -w 'files.docking.org -> HTTP %{http_code} in %{time_total}s\n' \
    https://files.docking.org/zinc22/ || echo "cannot reach files.docking.org"

# Machine-readable summary, read by the GUI's "Load server report" button
# (resources.parse_server_report). Keep the DF_ prefix and KEY=VALUE shape.
section "DOCKINGFLOW SUMMARY"
echo "DF_HOSTNAME=$(hostname)"
echo "DF_LOGICAL_CPUS=$(nproc 2>/dev/null || sysctl -n hw.logicalcpu)"
if have lscpu; then
    echo "DF_PHYSICAL_CORES=$(lscpu -p=CORE,SOCKET 2>/dev/null | grep -v '^#' | sort -u | wc -l)"
else
    echo "DF_PHYSICAL_CORES=$(sysctl -n hw.physicalcpu 2>/dev/null)"
fi
if [ -r /proc/meminfo ]; then
    echo "DF_MEM_TOTAL_KB=$(awk '/^MemTotal:/ {print $2}' /proc/meminfo)"
    echo "DF_MEM_AVAILABLE_KB=$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo)"
else
    echo "DF_MEM_TOTAL_KB=$(( $(sysctl -n hw.memsize) / 1024 ))"
fi
echo "DF_LOAD1=$(cut -d' ' -f1 /proc/loadavg 2>/dev/null || sysctl -n vm.loadavg | awk '{print $2}')"
echo "DF_DISK_FREE_KB=$(df -Pk "$WORKDIR" 2>/dev/null | awk 'NR==2 {print $4}')"
echo "DF_INODES_FREE=$(df -Pi "$WORKDIR" 2>/dev/null | awk 'NR==2 {print $4}')"

echo
echo "Done."
