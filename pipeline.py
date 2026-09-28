#!/usr/bin/env python3
"""DockingFlow pipeline orchestrator.

Runs a resumable, per-tranche virtual-screening pipeline against a ZINC
ligand library: download a tranche's compressed ligands, unpack them, dock
them against one or more (receptor, grid box) targets with a parallel
AutoDock Vina backend (e.g. VinaLC), then rank and filter the results.

Every tranche is a self-contained working directory under
``<workdir>/tranches/<label>/`` whose progress is tracked in a
``status.txt`` file, advancing through::

    INIT -> DOWNLOADED -> UNPACKED -> DOCKED -> DONE

(or ``FAILED_<STAGE>`` if that stage errors). Re-running the pipeline over
the same workdir skips any stage a tranche has already completed, so an
interrupted run (network hiccup, killed job, etc.) can simply be restarted.

See README.md for the full stage-by-stage description, configuration file
formats, and CLI reference. This module is also imported directly by
``gui.py``, which drives the same stage functions from a background thread
instead of ``main()``'s straight-through CLI loop.
"""
from __future__ import annotations

import argparse
import gzip
import shutil
import subprocess
import tarfile
import time
import zlib
from datetime import datetime
from pathlib import Path
from typing import Iterator

import io_parse
from io_parse import (
    DockingTarget,
    Tranche,
    VinaLCOptions,
    format_summary,
    load_docking_targets,
    load_setup,
    load_tranches_tsv,
    load_vinalc_options,
    validate_inputs,
)

# ----------------------------
# Small utilities
# ----------------------------

# Statuses a tranche passes through in order. FAILED_* statuses are terminal
# and deliberately excluded so a failed tranche is never treated as "past" a
# stage it never completed.
STATUS_ORDER = ["INIT", "DOWNLOADED", "UNPACKED", "DOCKED", "DONE"]


def status_at_least(status: str, target: str) -> bool:
    """Return whether `status` has reached or passed `target` in STATUS_ORDER.

    Any status not in STATUS_ORDER (in practice, a `FAILED_*` status) is
    treated as *not* having reached anything, so a failed tranche is always
    re-attempted from its failed stage rather than silently skipped.
    """
    try:
        return STATUS_ORDER.index(status) >= STATUS_ORDER.index(target)
    except ValueError:
        return False


tranche_label = io_parse.tranche_label


def write_text(path: Path, text: str) -> None:
    """Write `text` to `path`, creating parent directories as needed."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def read_status(tdir: Path) -> str:
    """Read a tranche's current status, defaulting to "INIT" if unset."""
    p = tdir / "status.txt"
    if not p.exists():
        return "INIT"
    return p.read_text(encoding="utf-8").strip() or "INIT"


def write_status(tdir: Path, status: str) -> None:
    """Persist a tranche's status. This is the single source of truth `read_status` reads back."""
    write_text(tdir / "status.txt", status + "\n")


def run_cmd(cmd: list[str], cwd: Path, stdout_path: Path, stderr_path: Path) -> tuple[int, float]:
    """Run `cmd` as a subprocess, redirecting stdout/stderr to files on disk.

    Used for both the download stage (curl script) and the docking stage
    (vinalc/mpirun invocation), so every external process this pipeline
    shells out to gets its own timed, on-disk log pair for post-mortem
    debugging.

    Args:
        cmd: Argv list passed straight to `subprocess.Popen` (no shell).
        cwd: Working directory the command is run from.
        stdout_path: File to write the command's stdout to.
        stderr_path: File to write the command's stderr to.

    Returns:
        A `(returncode, elapsed_seconds)` tuple.
    """
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    stderr_path.parent.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    with open(stdout_path, "w", encoding="utf-8") as out, open(stderr_path, "w", encoding="utf-8") as err:
        proc = subprocess.Popen(cmd, cwd=str(cwd), stdout=out, stderr=err, text=True)
        rc = proc.wait()
    return rc, time.time() - t0


# Downloaded ligand files the unpack stage knows how to read: ZINC22 3D
# archives (.pdbqt.tgz), ZINC20-style gzipped pdbqt, and plain pdbqt.
ARCHIVE_SUFFIXES = (".pdbqt.tgz", ".pdbqt.tar.gz", ".pdbqt.gz", ".pdbqt")


def find_ligand_archives(download_dir: Path) -> list[Path]:
    """Every downloaded ligand file under `download_dir`, recursively, in sorted order.

    ZINC downloads land in nested tranche subdirectories (e.g.
    `H04/H04M000/a/H04M000-N-aaaaaa.pdbqt.tgz`), so this must recurse rather
    than glob the top level only.
    """
    return sorted(p for p in download_dir.rglob("*") if p.is_file() and p.name.endswith(ARCHIVE_SUFFIXES))


def count_ligand_archives(download_dir: Path) -> int:
    return len(find_ligand_archives(download_dir))


# ----------------------------
# Stage 0: workspace creation
# ----------------------------

def materialize_tranche(workdir: Path, tranche: Tranche, setup_path: Path) -> Path:
    """Create (or reopen) a tranche's workspace directory and snapshot its inputs.

    Produces `<workdir>/tranches/<label>/` containing:
      - `inputs/curl_script.curl`, `inputs/setup_used.txt`: copies of the
        curl script and setup.txt used for this tranche, refreshed on every
        call (in case the mapping file's curl script changed) — kept as an
        audit trail so a run can be reproduced or debugged after the fact.
      - `status.txt`: created once, as "INIT", and never overwritten here.
      - `meta.txt`: created once, recording when the tranche was first
        materialized and its tranche identity.

    Safe to call repeatedly (e.g. on every pipeline re-run): it never
    overwrites `status.txt` or `meta.txt` once they exist, so a tranche's
    progress and creation time survive across runs.

    Args:
        workdir: The run's top-level output directory.
        tranche: The tranche to materialize.
        setup_path: Path to the setup.txt in effect for this run, copied
            into the tranche's `inputs/` for the audit trail.

    Returns:
        The tranche's workspace directory.
    """
    tdir = workdir / "tranches" / tranche_label(tranche)
    inputs_dir = tdir / "inputs"
    inputs_dir.mkdir(parents=True, exist_ok=True)

    # Snapshot artifacts for audit trail (always refresh; curl scripts may update)
    shutil.copy2(tranche.curl_script, inputs_dir / "curl_script.curl")
    shutil.copy2(setup_path, inputs_dir / "setup_used.txt")

    # Status/meta created once
    status_path = tdir / "status.txt"
    if not status_path.exists():
        write_status(tdir, "INIT")

    meta_path = tdir / "meta.txt"
    if not meta_path.exists():
        write_text(
            meta_path,
            f"created_at={datetime.now().isoformat(timespec='seconds')}\n"
            f"tranche={tranche_label(tranche)}\n"
            f"log_p={tranche.log_p}\n"
            f"size={io_parse.tranche_size(tranche)}\n"
            f"curl_script_original={tranche.curl_script}\n"
        )

    return tdir


# ----------------------------
# Stage 1: download
# ----------------------------

def download_tranche(tdir: Path) -> None:
    """Run a tranche's snapshotted curl script and verify it produced ligands.

    Skips entirely if the tranche has already reached "DOWNLOADED" or later
    (idempotent/resumable). If the script produced no ligand files at all,
    writes "FAILED_DOWNLOAD" and raises, so a caller looping over tranches
    can decide whether to abort (CLI) or record the error and continue with
    the next tranche (GUI).

    A nonzero exit code with *some* files downloaded is treated as a partial
    success: ZINC regularly 404s on individual archives, and failing the
    whole tranche on those would make it impossible to ever finish. The
    failed commands (as reported by the `zinc_split.py` script wrapper) are
    written to `<tdir>/download_failures.txt` and a warning is printed; to
    retry them, delete the tranche directory and re-run.

    Args:
        tdir: The tranche's workspace directory, as returned by
            `materialize_tranche`. Must contain
            `inputs/curl_script.curl`.

    Raises:
        FileNotFoundError: If the curl script snapshot is missing.
        RuntimeError: If the curl script downloads no ligand files.
    """
    status = read_status(tdir)
    if status_at_least(status, "DOWNLOADED"):
        print(f"[skip] {tdir.name}: status={status}")
        return

    print(f"[download] {tdir.name}")

    download_dir = tdir / "download"
    logs_dir = tdir / "logs"
    download_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    curl_script = (tdir / "inputs" / "curl_script.curl").resolve()  # absolute: it runs from download_dir
    if not curl_script.exists():
        write_status(tdir, "FAILED_DOWNLOAD")
        raise FileNotFoundError(f"Missing curl script snapshot: {curl_script}")

    stdout_log = logs_dir / "download.stdout.log"
    stderr_log = logs_dir / "download.stderr.log"

    # Run curl script in download_dir so --create-dirs writes underneath it.
    rc, dt = run_cmd(["bash", str(curl_script)], cwd=download_dir, stdout_path=stdout_log, stderr_path=stderr_log)

    n_files = count_ligand_archives(download_dir)
    write_text(tdir / "download_timing.txt", f"rc={rc}\nseconds={dt:.2f}\nfiles={n_files}\n")

    if n_files == 0:
        write_status(tdir, "FAILED_DOWNLOAD")
        raise RuntimeError(f"No ligand files downloaded for {tdir.name} (rc={rc}). See {stderr_log}")

    if rc != 0:
        failed = [ln for ln in stderr_log.read_text(errors="replace").splitlines() if ln.startswith("FAILED")]
        failures_path = tdir / "download_failures.txt"
        write_text(failures_path, "".join(f"{ln}\n" for ln in failed))
        print(
            f"[warn] {tdir.name}: download script exited rc={rc} with {len(failed)} failed command(s); "
            f"continuing with {n_files} downloaded file(s). See {failures_path}"
        )

    write_status(tdir, "DOWNLOADED")


# ----------------------------
# Stage 2: unpack
# ----------------------------

# Molecules per generated ligand file. VinaLC reads each file listed in its
# ligList sequentially, so this only bounds file size (a ZINC22 tranche can
# hold hundreds of thousands of molecules), not parallelism.
LIGANDS_PER_FILE = 10000

NAME_REMARK_PREFIX = "REMARK  Name = "


def iter_pdbqt_texts(archive: Path) -> Iterator[tuple[str, str]]:
    """Yield `(source_name, pdbqt_text)` for every pdbqt file in a downloaded file.

    Handles ZINC22 `.pdbqt.tgz` archives (one small pdbqt per molecule,
    read straight out of the archive without extracting to disk),
    ZINC20-style `.pdbqt.gz`, and plain `.pdbqt`.
    """
    if archive.name.endswith((".pdbqt.tgz", ".pdbqt.tar.gz")):
        with tarfile.open(archive, "r:gz") as tar:
            for member in tar:
                if not member.isfile() or not member.name.endswith(".pdbqt"):
                    continue
                f = tar.extractfile(member)
                if f is not None:
                    yield member.name, f.read().decode("utf-8", errors="replace")
    elif archive.name.endswith(".pdbqt.gz"):
        with gzip.open(archive, "rt", encoding="utf-8", errors="replace") as f:
            yield archive.name, f.read()
    else:
        yield archive.name, archive.read_text(encoding="utf-8", errors="replace")


def ligand_name(lines: list[str]) -> str | None:
    """The molecule's `REMARK  Name = <id>` value (ZINC writes one per molecule), if any."""
    for line in lines:
        if line.startswith("REMARK") and "Name =" in line:
            return line.split("=", 1)[1].strip() or None
    return None


def split_molecules(source_name: str, text: str) -> Iterator[tuple[str, list[str]]]:
    """Yield `(ligand_id, lines)` for each molecule in one pdbqt file's text.

    A file with `MODEL`/`ENDMDL` blocks holds one molecule per block; a file
    without them is a single molecule. The returned lines never include
    `MODEL`/`ENDMDL` themselves. The id is the molecule's
    `REMARK  Name = ...` value, falling back to the source file name.
    """
    stem = source_name.rsplit("/", 1)[-1]
    for suffix in (".gz", ".pdbqt"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]

    lines = text.splitlines()
    if any(ln.startswith("MODEL") for ln in lines):
        blocks: list[list[str]] = []
        current: list[str] | None = None
        for ln in lines:
            if ln.startswith("MODEL"):
                current = []
            elif ln.startswith("ENDMDL"):
                if current is not None:
                    blocks.append(current)
                current = None
            elif current is not None:
                current.append(ln)
    else:
        blocks = [lines]

    for i, block in enumerate(blocks, start=1):
        if not any(ln.startswith(("ATOM", "HETATM")) for ln in block):
            continue
        name = ligand_name(block) or (stem if len(blocks) == 1 else f"{stem}_{i}")
        yield name, block


def unpack_tranche(tdir: Path) -> Path:
    """Repackage a tranche's downloaded ligands into VinaLC's multi-molecule input format.

    VinaLC treats every `MODEL`...`ENDMDL` block of each file in its ligList
    as one docking job, and ignores anything outside those blocks — so a
    plain single-molecule pdbqt (which is what ZINC22 archives contain)
    would silently produce zero jobs. This reads every molecule out of
    every downloaded file and writes them, each wrapped in `MODEL`/`ENDMDL`,
    into `<tdir>/ligands/ligands_00001.pdbqt` etc. (`LIGANDS_PER_FILE` per
    file). It also writes:

      - `<tdir>/ligand_list.txt`: those files' paths, one per line, in
        order — VinaLC's `--ligList`.
      - `<tdir>/ligand_index.tsv`: `index  ligand  source`. VinaLC labels
        results only as `LIGAND <n>`, counting molecules across the ligList
        files in order, so this maps `n` back to the ZINC id.

    A downloaded file that can't be read (e.g. a truncated download) is
    skipped and listed in `<tdir>/unpack_warnings.txt`.

    Skips entirely if already "UNPACKED" or later. Requires the tranche to
    already be "DOWNLOADED".

    Args:
        tdir: The tranche's workspace directory.

    Returns:
        The tranche's `ligands/` directory (whether newly unpacked or
        already existing from a prior run).

    Raises:
        RuntimeError: If called before the download stage has completed, or
            if no molecules could be read from the downloaded files.
    """
    status = read_status(tdir)
    ligands_dir = tdir / "ligands"
    if status_at_least(status, "UNPACKED"):
        print(f"[skip] {tdir.name}: status={status}")
        return ligands_dir
    if not status_at_least(status, "DOWNLOADED"):
        raise RuntimeError(f"{tdir.name}: cannot unpack before download (status={status})")

    print(f"[unpack] {tdir.name}")

    archives = find_ligand_archives(tdir / "download")
    if not archives:
        write_status(tdir, "FAILED_UNPACK")
        raise RuntimeError(f"No ligand files found to unpack in {tdir / 'download'}")

    # Start clean so a previously interrupted unpack can't leave stale chunks behind.
    if ligands_dir.exists():
        shutil.rmtree(ligands_dir)
    ligands_dir.mkdir(parents=True)

    ligand_files: list[Path] = []
    warnings: list[str] = []
    n = 0
    out = None
    with open(tdir / "ligand_index.tsv", "w", encoding="utf-8") as index:
        index.write("index\tligand\tsource\n")
        try:
            for archive in archives:
                try:
                    for source_name, text in iter_pdbqt_texts(archive):
                        for name, lines in split_molecules(source_name, text):
                            if n % LIGANDS_PER_FILE == 0:
                                if out:
                                    out.close()
                                ligand_files.append(ligands_dir / f"ligands_{len(ligand_files) + 1:05d}.pdbqt")
                                out = open(ligand_files[-1], "w", encoding="utf-8")
                            n += 1
                            out.write(f"MODEL {n}\n")
                            if ligand_name(lines) is None:
                                # Carry the id into VinaLC's output poses too.
                                out.write(f"{NAME_REMARK_PREFIX}{name}\n")
                            out.write("\n".join(lines) + "\nENDMDL\n")
                            index.write(f"{n}\t{name}\t{archive.name}\n")
                except (tarfile.TarError, OSError, EOFError, zlib.error) as exc:
                    warnings.append(f"{archive}: {exc}")
        finally:
            if out:
                out.close()

    write_text(tdir / "unpack_warnings.txt", "".join(f"{w}\n" for w in warnings))
    if warnings:
        print(f"[warn] {tdir.name}: skipped {len(warnings)} unreadable file(s). See {tdir / 'unpack_warnings.txt'}")

    if n == 0:
        write_status(tdir, "FAILED_UNPACK")
        raise RuntimeError(f"{tdir.name}: no molecules found in {len(archives)} downloaded file(s)")

    print(f"[unpack] {tdir.name}: {n} molecule(s) from {len(archives)} file(s) into {len(ligand_files)} ligand file(s)")
    write_text(tdir / "ligand_list.txt", "".join(f"{p}\n" for p in ligand_files))
    write_status(tdir, "UNPACKED")
    return ligands_dir


# ----------------------------
# Stage 3: dock
# ----------------------------

# VinaLC writes its combined results to "<recList arg>_<ligList arg>.pdbqt.gz"
# (and ".log.gz") in its working directory, using the arguments verbatim —
# so it's always run from the docking dir with these bare file names.
VINALC_REC_LIST = "recList.txt"
VINALC_GEO_LIST = "geoList.txt"
VINALC_LIG_LIST = "ligList.txt"
VINALC_POSES = f"{VINALC_REC_LIST}_{VINALC_LIG_LIST}.pdbqt.gz"


def write_vinalc_inputs(docking_dir: Path, targets: list[DockingTarget], ligand_list_path: Path) -> None:
    """Write the three list files VinaLC reads, with absolute paths, into `docking_dir`."""
    write_text(docking_dir / VINALC_REC_LIST, "".join(f"{t.receptor}\n" for t in targets))
    write_text(
        docking_dir / VINALC_GEO_LIST,
        "".join(
            f"{b.center_x} {b.center_y} {b.center_z} {b.size_x} {b.size_y} {b.size_z}\n"
            for b in (t.box for t in targets)
        ),
    )
    shutil.copyfile(ligand_list_path, docking_dir / VINALC_LIG_LIST)


def resolve_bin(name: str) -> str:
    """Make a relative binary path (e.g. `./bin/vinalc`) absolute, since docking runs from another dir.

    Bare names (`vinalc`, `mpirun`) are left for PATH lookup.
    """
    return str(Path(name).expanduser().resolve()) if "/" in name else name


def build_vinalc_command(options: VinaLCOptions, vinalc_bin: str, mpirun_bin: str | None) -> list[str]:
    """The full docking command line, run from the tranche's docking dir."""
    cmd = [
        resolve_bin(vinalc_bin),
        "--recList", VINALC_REC_LIST,
        "--ligList", VINALC_LIG_LIST,
        "--geoList", VINALC_GEO_LIST,
        *options.cli_args(),
    ]
    if mpirun_bin:
        cmd = [resolve_bin(mpirun_bin), "-np", str(options.mpi_ranks), *cmd]
    return cmd


def dock_tranche(
    tdir: Path,
    targets: list[DockingTarget],
    *,
    options: VinaLCOptions,
    vinalc_bin: str,
    mpirun_bin: str | None,
) -> None:
    """Dock a tranche's unpacked ligands against every configured target with one VinaLC run.

    VinaLC itself loops over every (receptor, grid box) in its recList/
    geoList and every molecule in its ligList, farming each pair out to its
    MPI worker ranks, so the whole tranche is a single
    `mpirun -np <ranks> vinalc ...` call run from `<tdir>/docking/`. Its
    combined poses land in `<tdir>/docking/recList.txt_ligList.txt.pdbqt.gz`
    (with a matching `.log.gz`).

    Skips entirely if already "DOCKED" or later. Requires "UNPACKED".

    Args:
        tdir: The tranche's workspace directory.
        targets: One or more `(receptor, grid box)` pairs, as returned by
            `io_parse.load_docking_targets`.
        options: Parsed docking settings (`io_parse.load_vinalc_options`),
            including the MPI rank count.
        vinalc_bin: Path/name of the VinaLC binary.
        mpirun_bin: MPI launcher to wrap `vinalc_bin` with, or `None` to
            invoke `vinalc_bin` directly (only useful for the offline test
            mock; real VinaLC exits unless it has at least 2 MPI ranks).

    Raises:
        RuntimeError: If called before unpacking, if `ligand_list.txt` is
            missing/empty, or if VinaLC exits nonzero or writes no output.
    """
    status = read_status(tdir)
    if status_at_least(status, "DOCKED"):
        print(f"[skip] {tdir.name}: status={status}")
        return
    if not status_at_least(status, "UNPACKED"):
        raise RuntimeError(f"{tdir.name}: cannot dock before unpack (status={status})")

    print(f"[dock] {tdir.name}")

    ligand_list_path = tdir / "ligand_list.txt"
    if not ligand_list_path.exists() or not ligand_list_path.read_text().strip():
        write_status(tdir, "FAILED_DOCK")
        raise RuntimeError(f"{tdir.name}: no ligand_list.txt to dock")

    docking_dir = tdir / "docking"
    logs_dir = tdir / "logs"
    docking_dir.mkdir(parents=True, exist_ok=True)
    poses_path = docking_dir / VINALC_POSES
    poses_path.unlink(missing_ok=True)  # never rank a stale result from an earlier failed attempt

    write_vinalc_inputs(docking_dir, targets, ligand_list_path)
    cmd = build_vinalc_command(options, vinalc_bin, mpirun_bin)
    write_text(docking_dir / "command.txt", " ".join(cmd) + "\n")

    stdout_log = logs_dir / "dock.stdout.log"
    stderr_log = logs_dir / "dock.stderr.log"
    try:
        rc, dt = run_cmd(cmd, cwd=docking_dir, stdout_path=stdout_log, stderr_path=stderr_log)
    except FileNotFoundError:
        write_status(tdir, "FAILED_DOCK")
        raise RuntimeError(f"{tdir.name}: docking binary not found: {cmd[0]} (put it on PATH or give its full path)")
    write_text(docking_dir / "dock_timing.txt", f"rc={rc}\nseconds={dt:.2f}\n")

    if rc != 0:
        write_status(tdir, "FAILED_DOCK")
        raise RuntimeError(f"Docking failed for {tdir.name} (rc={rc}). See {stderr_log}")
    if not poses_path.exists():
        write_status(tdir, "FAILED_DOCK")
        raise RuntimeError(f"Docking for {tdir.name} exited 0 but wrote no {poses_path.name}. See {stdout_log}")

    write_status(tdir, "DOCKED")


# ----------------------------
# Stage 4: rank + filter
# ----------------------------

def parse_vinalc_poses(poses_path: Path) -> Iterator[tuple[str, int, float]]:
    """Yield `(receptor, ligand_number, best_affinity)` per docked pair in a VinaLC output file.

    VinaLC appends one record per (receptor, ligand) job, in completion
    order (not input order)::

        REMARK RECEPTOR /abs/path/protein.pdbqt
        REMARK LIGAND LIGAND 12
        MODEL 1
        REMARK VINA RESULT:      -7.4      0.000      0.000
        ...
        ENDMDL
        MODEL 2 ...

    Poses are best-first, so the first `REMARK VINA RESULT:` of a record is
    its best affinity. A job that found no pose inside the box has no
    `VINA RESULT` line and is skipped.
    """
    receptor: str | None = None
    number: int | None = None
    best: float | None = None
    with gzip.open(poses_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("REMARK RECEPTOR"):
                if receptor is not None and number is not None and best is not None:
                    yield receptor, number, best
                receptor, number, best = line[len("REMARK RECEPTOR"):].strip(), None, None
            elif line.startswith("REMARK LIGAND"):
                number = int(line.split()[-1])
            elif best is None and line.startswith("REMARK VINA RESULT:"):
                best = float(line.split()[3])
    if receptor is not None and number is not None and best is not None:
        yield receptor, number, best


def read_ligand_index(tdir: Path) -> dict[int, str]:
    """Map VinaLC's `LIGAND <n>` numbers back to ligand ids, from `ligand_index.tsv`."""
    index: dict[int, str] = {}
    with open(tdir / "ligand_index.tsv", encoding="utf-8") as f:
        next(f)  # header
        for line in f:
            n, name, _source = line.rstrip("\n").split("\t")
            index[int(n)] = name
    return index


def rank_and_filter_tranche(tdir: Path, filter_percent: float) -> Path:
    """Rank a tranche's docked ligands by best binding affinity and keep the top slice.

    Reads VinaLC's combined output, maps each `LIGAND <n>` back to its ZINC
    id via `ligand_index.tsv`, keeps each ligand's best affinity across all
    receptors (noting which receptor it came from), sorts ascending (more
    negative kcal/mol = stronger predicted binding = better), and writes:

      - `<tdir>/results/ranked_all.tsv`: every docked ligand, best-first.
      - `<tdir>/results/top_hits.tsv`: the top `filter_percent`% slice
        (at least 1 ligand, even if `filter_percent` would round down to 0).
      - `<tdir>/results/summary.txt`: how many molecules were docked
        successfully out of how many were submitted.

    Requires the tranche to already be "DOCKED"; sets it to "DONE" on success.

    Args:
        tdir: The tranche's workspace directory.
        filter_percent: Percentage (0, 100] of ranked ligands to keep in
            `top_hits.tsv`.

    Returns:
        Path to `top_hits.tsv`.

    Raises:
        RuntimeError: If called before docking, or if no docking results
            with a parseable affinity are found.
    """
    status = read_status(tdir)
    if not status_at_least(status, "DOCKED"):
        raise RuntimeError(f"{tdir.name}: cannot rank before docking (status={status})")

    index = read_ligand_index(tdir)
    best: dict[str, tuple[float, str]] = {}
    for receptor, number, affinity in parse_vinalc_poses(tdir / "docking" / VINALC_POSES):
        lig = index.get(number, f"LIGAND_{number}")
        if lig not in best or affinity < best[lig][0]:
            best[lig] = (affinity, Path(receptor).stem)

    if not best:
        write_status(tdir, "FAILED_RANK")
        raise RuntimeError(f"{tdir.name}: no docking results found to rank")

    ranked = sorted(best.items(), key=lambda kv: kv[1][0])
    header = "ligand\taffinity_kcal_mol\treceptor\n"
    results_dir = tdir / "results"
    write_text(
        results_dir / "ranked_all.tsv",
        header + "".join(f"{lig}\t{aff:.3f}\t{rec}\n" for lig, (aff, rec) in ranked),
    )

    n_keep = max(1, round(len(ranked) * filter_percent / 100))
    top_path = results_dir / "top_hits.tsv"
    write_text(top_path, header + "".join(f"{lig}\t{aff:.3f}\t{rec}\n" for lig, (aff, rec) in ranked[:n_keep]))

    submitted = len(set(index.values()))
    write_text(results_dir / "summary.txt", f"submitted={submitted}\ndocked={len(ranked)}\ntop_hits={n_keep}\n")
    if len(ranked) < submitted:
        print(f"[warn] {tdir.name}: only {len(ranked)} of {submitted} ligands produced a pose")

    write_status(tdir, "DONE")
    return top_path


def combine_results(workdir: Path, tranche_dirs: list[Path]) -> Path:
    """Merge every tranche's top_hits.tsv into one workdir-level ranking."""
    combined: list[tuple[str, float, str, str]] = []
    for tdir in tranche_dirs:
        top_path = tdir / "results" / "top_hits.tsv"
        if not top_path.exists():
            continue
        for line in top_path.read_text().splitlines()[1:]:
            if not line.strip():
                continue
            lig, aff, receptor = line.split("\t")
            combined.append((lig, float(aff), receptor, tdir.name))

    combined.sort(key=lambda row: row[1])
    out_path = workdir / "results" / "top_hits_combined.tsv"
    write_text(
        out_path,
        "ligand\taffinity_kcal_mol\treceptor\ttranche\n"
        + "".join(f"{lig}\t{aff:.3f}\t{rec}\t{tranche}\n" for lig, aff, rec, tranche in combined),
    )
    return out_path


# ----------------------------
# Cleaning
# ----------------------------

def clean_workdir(workdir: Path) -> None:
    """
    Delete all generated tranche directories under workdir/tranches/.
    Intended for test resets. Keeps the workdir itself.
    """
    tranches_dir = workdir / "tranches"
    if not tranches_dir.exists():
        print("Nothing to clean (no tranches directory).")
        return

    # Safety guard: never delete root or home by accident
    if tranches_dir == Path("/") or tranches_dir == Path.home():
        raise RuntimeError(f"Refusing to delete dangerous path: {tranches_dir}")

    shutil.rmtree(tranches_dir)
    print(f"Deleted generated tranche directories: {tranches_dir}")


def nuke_workdir(workdir: Path) -> None:
    """Delete the entire workdir. DANGEROUS: intended for full resets."""
    if not workdir.exists():
        print(f"Nothing to nuke (workdir does not exist): {workdir}")
        return

    if workdir == Path("/") or workdir == Path.home():
        raise RuntimeError(f"Refusing to delete dangerous path: {workdir}")

    shutil.rmtree(workdir)
    print(f"Deleted workdir: {workdir}")

# ----------------------------
# Main
# ----------------------------

def main() -> int:
    """CLI entry point: parse arguments and run the pipeline stages in order.

    Unlike the GUI's runner (`gui.PipelineAPI._run`), this stops the whole
    run on the first tranche that raises — there's no per-tranche error
    isolation here, on the theory that an unattended CLI run should fail
    loudly rather than silently produce partial results. Use
    `--only-stage0` / `--only-download` / `--only-unpack` to stop early for
    inspection, or `--clean` / `--nuke` to reset a workdir between attempts.
    """
    parser = argparse.ArgumentParser(
        description="DockingFlow: download ZINC tranches, dock them with VinaLC, rank the hits"
    )
    parser.add_argument("--map", required=True, help="Path to tranches.tsv (tab-separated)")
    parser.add_argument("--setup", required=True, help="Path to setup.txt (KEY=VALUE)")
    parser.add_argument("--workdir", required=True, help="Output directory for this run")
    parser.add_argument("--clean", action="store_true", help="Delete generated outputs in workdir and reset statuses")
    parser.add_argument("--nuke", action="store_true", help="Delete the entire workdir (DANGEROUS)")
    parser.add_argument("--only-stage0", action="store_true", help="Only create tranche workspaces (no downloads)")
    parser.add_argument("--only-download", action="store_true", help="Stop after the download stage")
    parser.add_argument("--only-unpack", action="store_true", help="Stop after unpacking downloaded ligands")
    parser.add_argument("--vinalc-bin", default="vinalc", help="Docking binary to invoke (default: vinalc)")
    parser.add_argument("--mpirun-bin", default="mpirun", help="MPI launcher to wrap the docking binary with")
    parser.add_argument(
        "--no-mpirun", action="store_true", help="Invoke the docking binary directly, without an MPI launcher"
    )
    args = parser.parse_args()

    workdir = Path(args.workdir).expanduser().resolve()

    # Clean/nuke modes should not require valid inputs.
    if args.nuke:
        nuke_workdir(workdir)
        return 0

    workdir.mkdir(parents=True, exist_ok=True)

    if args.clean:
        clean_workdir(workdir)
        return 0

    setup_path = Path(args.setup).expanduser().resolve()
    map_path = Path(args.map).expanduser().resolve()

    setup = load_setup(setup_path)
    tranches = load_tranches_tsv(map_path)

    validate_inputs(setup, tranches, workdir)

    print(format_summary(setup, tranches, workdir))

    (workdir / "tranches").mkdir(exist_ok=True)

    # Stage 0: create tranche directories + snapshot inputs
    tranche_dirs: list[Path] = []
    for t in tranches:
        tdir = materialize_tranche(workdir, t, setup_path)
        tranche_dirs.append(tdir)

    print(f"\nStage 0 complete: prepared {len(tranche_dirs)} tranche workspaces.")

    if args.only_stage0:
        return 0

    # Stage 1: download each tranche sequentially (safe default)
    for tdir in tranche_dirs:
        download_tranche(tdir)

    print("\nDownload stage complete for all tranches.")
    if args.only_download:
        return 0

    # Stage 2: unpack downloaded ligands
    for tdir in tranche_dirs:
        unpack_tranche(tdir)

    print("\nUnpack stage complete for all tranches.")
    if args.only_unpack:
        return 0

    # Stage 3: dock each tranche's ligands against every (receptor, grid box) target
    targets = load_docking_targets(setup)
    options = load_vinalc_options(setup)
    mpirun_bin = None if args.no_mpirun else args.mpirun_bin

    for tdir in tranche_dirs:
        dock_tranche(
            tdir,
            targets,
            options=options,
            vinalc_bin=args.vinalc_bin,
            mpirun_bin=mpirun_bin,
        )

    print("\nDocking stage complete for all tranches.")

    # Stage 4: rank + filter results, then combine across tranches
    filter_percent = float(setup["filter_percent"])
    for tdir in tranche_dirs:
        rank_and_filter_tranche(tdir, filter_percent)

    combined_path = combine_results(workdir, tranche_dirs)
    print(f"\nPipeline complete for all tranches. Combined top hits: {combined_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
