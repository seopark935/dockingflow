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
import time
from datetime import datetime
from pathlib import Path

from io_parse import (
    DockingTarget,
    Tranche,
    format_summary,
    load_docking_targets,
    load_setup,
    load_tranches_tsv,
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


def tranche_label(t: Tranche) -> str:
    """Filesystem-safe, deterministic tranche label."""
    logp = f"{t.log_p:.2f}".rstrip("0").rstrip(".")
    return f"LP{logp}_MW{t.molecular_weight}"


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


def count_pdbqt_gz(download_dir: Path) -> int:
    """Count `*.pdbqt.gz` files anywhere under `download_dir`, recursively.

    ZINC downloads land one compound per subdirectory (e.g.
    `ZINC000000000001/mol1.pdbqt.gz`), so this must recurse rather than glob
    the top level only.
    """
    return sum(1 for _ in download_dir.rglob("*.pdbqt.gz"))


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
        materialized and its (log_p, molecular_weight) identity.

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
            f"log_p={tranche.log_p}\n"
            f"molecular_weight={tranche.molecular_weight}\n"
            f"curl_script_original={tranche.curl_script}\n"
        )

    return tdir


# ----------------------------
# Stage 1: download
# ----------------------------

def download_tranche(tdir: Path) -> None:
    """Run a tranche's snapshotted curl script and verify it produced ligands.

    Skips entirely if the tranche has already reached "DOWNLOADED" or later
    (idempotent/resumable). On failure — nonzero exit code, or an exit code
    of 0 but zero `*.pdbqt.gz` files produced — writes "FAILED_DOWNLOAD" and
    raises, so a caller looping over tranches can decide whether to abort
    (CLI) or record the error and continue with the next tranche (GUI).

    Args:
        tdir: The tranche's workspace directory, as returned by
            `materialize_tranche`. Must contain
            `inputs/curl_script.curl`.

    Raises:
        FileNotFoundError: If the curl script snapshot is missing.
        RuntimeError: If the curl script exits nonzero, or exits 0 but
            downloads no ligand files.
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

    curl_script = tdir / "inputs" / "curl_script.curl"
    if not curl_script.exists():
        write_status(tdir, "FAILED_DOWNLOAD")
        raise FileNotFoundError(f"Missing curl script snapshot: {curl_script}")

    stdout_log = logs_dir / "download.stdout.log"
    stderr_log = logs_dir / "download.stderr.log"

    # Run curl script in download_dir so --create-dirs writes underneath it.
    rc, dt = run_cmd(["bash", str(curl_script)], cwd=download_dir, stdout_path=stdout_log, stderr_path=stderr_log)

    n_files = count_pdbqt_gz(download_dir)
    write_text(tdir / "download_timing.txt", f"rc={rc}\nseconds={dt:.2f}\nfiles={n_files}\n")

    if rc != 0:
        write_status(tdir, "FAILED_DOWNLOAD")
        raise RuntimeError(f"Download failed for {tdir.name} (rc={rc}). See {stderr_log}")

    if n_files == 0:
        write_status(tdir, "FAILED_DOWNLOAD")
        raise RuntimeError(f"No *.pdbqt.gz files downloaded for {tdir.name}. See logs in {logs_dir}")

    write_status(tdir, "DOWNLOADED")


# ----------------------------
# Stage 2: unpack
# ----------------------------

def unpack_tranche(tdir: Path) -> Path:
    """Decompress a tranche's downloaded `*.pdbqt.gz` files into a flat ligand list.

    Every ligand file, regardless of which per-compound subdirectory it was
    downloaded into, is decompressed into a single flat `<tdir>/ligands/`
    directory, and its path is written (one per line) to
    `<tdir>/ligand_list.txt` — the file the docking stage points its docking
    binary's `ligand_list` config option at.

    Skips entirely if already "UNPACKED" or later. Requires the tranche to
    already be "DOWNLOADED".

    Args:
        tdir: The tranche's workspace directory.

    Returns:
        The tranche's `ligands/` directory (whether newly unpacked or
        already existing from a prior run).

    Raises:
        RuntimeError: If called before the download stage has completed, or
            if no `*.pdbqt.gz` files are found to unpack.
    """
    status = read_status(tdir)
    ligands_dir = tdir / "ligands"
    if status_at_least(status, "UNPACKED"):
        print(f"[skip] {tdir.name}: status={status}")
        return ligands_dir
    if not status_at_least(status, "DOWNLOADED"):
        raise RuntimeError(f"{tdir.name}: cannot unpack before download (status={status})")

    print(f"[unpack] {tdir.name}")

    download_dir = tdir / "download"
    ligands_dir.mkdir(parents=True, exist_ok=True)

    gz_files = sorted(download_dir.rglob("*.pdbqt.gz"))
    if not gz_files:
        write_status(tdir, "FAILED_UNPACK")
        raise RuntimeError(f"No *.pdbqt.gz files found to unpack in {download_dir}")

    ligand_paths: list[Path] = []
    for gz_path in gz_files:
        out_path = ligands_dir / gz_path.with_suffix("").name
        with gzip.open(gz_path, "rb") as src, open(out_path, "wb") as dst:
            shutil.copyfileobj(src, dst)
        ligand_paths.append(out_path)

    write_text(tdir / "ligand_list.txt", "\n".join(str(p) for p in ligand_paths) + "\n")
    write_status(tdir, "UNPACKED")
    return ligands_dir


# ----------------------------
# Stage 3: dock
# ----------------------------

def build_vinalc_config(
    out_dir: Path, receptor: Path, box, ligand_list_path: Path, energy_range: str, num_modes: int = 9
) -> Path:
    """Write a vinalc-style KEY = VALUE config file for one (receptor, grid box) target.

    This mirrors the standard AutoDock Vina config format (which VinaLC also
    accepts), so it should work unmodified against most Vina-family docking
    binaries. If your specific binary expects a different config dialect,
    this is the one function to adapt.

    Args:
        out_dir: Directory the config (and, once run, the docking outputs)
            will live in.
        receptor: Path to the receptor `.pdbqt`.
        box: A `GridBox`-shaped object with center_x/y/z and size_x/y/z.
        ligand_list_path: Path to the tranche's `ligand_list.txt`.
        energy_range: Value of setup.txt's `energy_range`, passed through
            verbatim as a string (Vina/VinaLC parses it itself).
        num_modes: Maximum number of binding poses to report per ligand.

    Returns:
        Path to the written config file.
    """
    cfg_path = out_dir / "vinalc.conf"
    write_text(
        cfg_path,
        f"receptor = {receptor}\n"
        f"ligand_list = {ligand_list_path}\n"
        f"center_x = {box.center_x}\n"
        f"center_y = {box.center_y}\n"
        f"center_z = {box.center_z}\n"
        f"size_x = {box.size_x}\n"
        f"size_y = {box.size_y}\n"
        f"size_z = {box.size_z}\n"
        f"energy_range = {energy_range}\n"
        f"num_modes = {num_modes}\n"
        f"out_dir = {out_dir}\n",
    )
    return cfg_path


def dock_tranche(
    tdir: Path,
    targets: list[DockingTarget],
    *,
    energy_range: str,
    cores: int,
    vinalc_bin: str,
    mpirun_bin: str | None,
) -> None:
    """Dock a tranche's unpacked ligands against every configured target.

    Runs one docking-binary invocation per `(receptor, grid box)` target in
    `targets`, each against the tranche's full `ligand_list.txt`. Each
    target gets its own subdirectory under `<tdir>/docking/target<i>_<receptor
    stem>/`, containing its config file, docked output `*_out.pdbqt` files,
    and timing info — so results from multiple targets never collide, and a
    later `rank_and_filter_tranche` call can be pointed at the merged best
    affinity across all of them.

    Skips entirely if already "DOCKED" or later. Requires "UNPACKED".

    Args:
        tdir: The tranche's workspace directory.
        targets: One or more `(receptor, grid box)` pairs, as returned by
            `io_parse.load_docking_targets`.
        energy_range: Passed straight through to `build_vinalc_config`.
        cores: Degree of parallelism; passed as `mpirun -np <cores>` when
            `mpirun_bin` is set.
        vinalc_bin: Path/name of the docking binary to invoke.
        mpirun_bin: MPI launcher to wrap `vinalc_bin` with, or `None` to
            invoke `vinalc_bin` directly (e.g. for a single-core run, or
            when testing against a mock binary that isn't MPI-aware).

    Raises:
        RuntimeError: If called before unpacking, if `ligand_list.txt` is
            missing/empty, or if the docking binary exits nonzero for any
            target.
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

    for i, target in enumerate(targets, start=1):
        target_dir = docking_dir / f"target{i}_{target.receptor.stem}"
        target_dir.mkdir(parents=True, exist_ok=True)

        cfg_path = build_vinalc_config(target_dir, target.receptor, target.box, ligand_list_path, energy_range)

        cmd = (
            [mpirun_bin, "-np", str(cores), vinalc_bin, "--config", str(cfg_path)]
            if mpirun_bin
            else [vinalc_bin, "--config", str(cfg_path)]
        )

        stdout_log = logs_dir / f"dock_target{i}.stdout.log"
        stderr_log = logs_dir / f"dock_target{i}.stderr.log"
        rc, dt = run_cmd(cmd, cwd=target_dir, stdout_path=stdout_log, stderr_path=stderr_log)
        write_text(target_dir / "dock_timing.txt", f"rc={rc}\nseconds={dt:.2f}\n")

        if rc != 0:
            write_status(tdir, "FAILED_DOCK")
            raise RuntimeError(f"Docking failed for {tdir.name} target {i} (rc={rc}). See {stderr_log}")

    write_status(tdir, "DOCKED")


# ----------------------------
# Stage 4: rank + filter
# ----------------------------

def parse_best_affinities(out_dir: Path) -> dict[str, float]:
    """Read docked *_out.pdbqt files and return {ligand_stem: best_affinity_kcal_mol}.

    Vina/VinaLC list poses best-first, so the first 'REMARK VINA RESULT:' line
    in each output file is that ligand's best-scoring pose.
    """
    best: dict[str, float] = {}
    for out_path in sorted(out_dir.glob("*_out.pdbqt")):
        for line in out_path.read_text().splitlines():
            if line.startswith("REMARK VINA RESULT:"):
                best[out_path.stem[: -len("_out")]] = float(line.split()[3])
                break
    return best


def rank_and_filter_tranche(tdir: Path, filter_percent: float) -> Path:
    """Rank a tranche's docked ligands by best binding affinity and keep the top slice.

    Merges best-affinity results across every target the tranche was docked
    against (if a ligand appears under multiple targets, the last target
    processed wins — in practice tranches are docked against a small, fixed
    set of targets, so this is rarely a meaningful ambiguity), sorts
    ascending by affinity (more negative kcal/mol = stronger predicted
    binding = better), and writes:

      - `<tdir>/results/ranked_all.tsv`: every docked ligand, best-first.
      - `<tdir>/results/top_hits.tsv`: the top `filter_percent`% slice
        (at least 1 ligand, even if `filter_percent` would round down to 0).

    Requires the tranche to already be "DOCKED"; sets it to "DONE" on success.

    Args:
        tdir: The tranche's workspace directory.
        filter_percent: Percentage (0, 100] of ranked ligands to keep in
            `top_hits.tsv`.

    Returns:
        Path to `top_hits.tsv`.

    Raises:
        RuntimeError: If called before docking, or if no docking results
            (`*_out.pdbqt` files with a parseable affinity) are found.
    """
    status = read_status(tdir)
    if not status_at_least(status, "DOCKED"):
        raise RuntimeError(f"{tdir.name}: cannot rank before docking (status={status})")

    affinities: dict[str, float] = {}
    for target_dir in sorted((tdir / "docking").glob("target*")):
        affinities.update(parse_best_affinities(target_dir))

    if not affinities:
        write_status(tdir, "FAILED_RANK")
        raise RuntimeError(f"{tdir.name}: no docking results found to rank")

    ranked = sorted(affinities.items(), key=lambda kv: kv[1])
    results_dir = tdir / "results"
    write_text(
        results_dir / "ranked_all.tsv",
        "ligand\taffinity_kcal_mol\n" + "\n".join(f"{lig}\t{aff:.3f}" for lig, aff in ranked) + "\n",
    )

    n_keep = max(1, round(len(ranked) * filter_percent / 100))
    top = ranked[:n_keep]
    top_path = results_dir / "top_hits.tsv"
    write_text(
        top_path,
        "ligand\taffinity_kcal_mol\n" + "\n".join(f"{lig}\t{aff:.3f}" for lig, aff in top) + "\n",
    )

    write_status(tdir, "DONE")
    return top_path


def combine_results(workdir: Path, tranche_dirs: list[Path]) -> Path:
    """Merge every tranche's top_hits.tsv into one workdir-level ranking."""
    combined: list[tuple[str, float, str]] = []
    for tdir in tranche_dirs:
        top_path = tdir / "results" / "top_hits.tsv"
        if not top_path.exists():
            continue
        for line in top_path.read_text().splitlines()[1:]:
            if not line.strip():
                continue
            lig, aff = line.split("\t")
            combined.append((lig, float(aff), tdir.name))

    combined.sort(key=lambda row: row[1])
    out_path = workdir / "results" / "top_hits_combined.tsv"
    write_text(
        out_path,
        "ligand\taffinity_kcal_mol\ttranche\n"
        + "\n".join(f"{lig}\t{aff:.3f}\t{tranche}" for lig, aff, tranche in combined)
        + "\n",
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
        description="Docking pipeline orchestrator (Stage 0 + Download stage)"
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
    cores = int(setup["cores"])
    energy_range = setup.get("energy_range", "3")
    mpirun_bin = None if args.no_mpirun else args.mpirun_bin

    for tdir in tranche_dirs:
        dock_tranche(
            tdir,
            targets,
            energy_range=energy_range,
            cores=cores,
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
