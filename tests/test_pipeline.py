#!/usr/bin/env python3
"""End-to-end tests for dockingflow, run entirely offline with a mock vinalc.

Run with: python3 -m unittest discover -s tests -v
(from the repo root, with the repo root on PYTHONPATH -- see run_tests.sh)
"""
from __future__ import annotations

import shutil
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

import io_parse
import pipeline

MOCK_VINALC = REPO_ROOT / "tests" / "fixtures" / "mock_vinalc.py"


class DockingFlowEndToEndTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="dockingflow_test_"))
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)
        self.workdir = self.tmp / "workdir"
        self.workdir.mkdir()

        self.setup_path = REPO_ROOT / "setup.txt"
        self.map_path = REPO_ROOT / "tranches.txt"

    def _load(self):
        setup = io_parse.load_setup(self.setup_path)
        tranches = io_parse.load_tranches_tsv(self.map_path)
        io_parse.validate_inputs(setup, tranches, self.workdir)
        return setup, tranches

    def test_validate_inputs_accepts_repo_fixtures(self):
        setup, tranches = self._load()
        self.assertEqual(len(tranches), 1)
        self.assertEqual(tranches[0].log_p, 5.0)
        self.assertEqual(tranches[0].molecular_weight, 200)

    def test_load_docking_targets_pairs_rec_and_geo(self):
        setup, _ = self._load()
        targets = io_parse.load_docking_targets(setup)
        self.assertEqual(len(targets), 1)
        self.assertTrue(targets[0].receptor.name.endswith("protein.pdbqt"))
        self.assertAlmostEqual(targets[0].box.center_x, 97.392)
        self.assertAlmostEqual(targets[0].box.size_z, 80.0)

    def test_load_docking_targets_rejects_mismatched_lengths(self):
        setup, _ = self._load()
        geo_path = self.tmp / "geoList_bad.txt"
        geo_path.write_text("1 2 3 4 5 6\n1 2 3 4 5 6\n")
        bad_setup = dict(setup)
        bad_setup["geoList"] = str(geo_path)
        with self.assertRaises(ValueError):
            io_parse.load_docking_targets(bad_setup)

    def test_full_pipeline_runs_end_to_end_with_mock_vinalc(self):
        setup, tranches = self._load()
        targets = io_parse.load_docking_targets(setup)

        tranche_dirs = []
        for t in tranches:
            tdir = pipeline.materialize_tranche(self.workdir, t, self.setup_path)
            self.assertEqual(pipeline.read_status(tdir), "INIT")
            tranche_dirs.append(tdir)

        for tdir in tranche_dirs:
            pipeline.download_tranche(tdir)
            self.assertEqual(pipeline.read_status(tdir), "DOWNLOADED")
            self.assertGreater(pipeline.count_pdbqt_gz(tdir / "download"), 0)

        for tdir in tranche_dirs:
            pipeline.unpack_tranche(tdir)
            self.assertEqual(pipeline.read_status(tdir), "UNPACKED")
            ligand_list = tdir / "ligand_list.txt"
            self.assertTrue(ligand_list.exists())
            self.assertEqual(len(ligand_list.read_text().splitlines()), 2)

        for tdir in tranche_dirs:
            pipeline.dock_tranche(
                tdir,
                targets,
                energy_range=setup.get("energy_range", "3"),
                cores=2,
                vinalc_bin=str(MOCK_VINALC),
                mpirun_bin=None,
            )
            self.assertEqual(pipeline.read_status(tdir), "DOCKED")

        filter_percent = float(setup["filter_percent"])
        for tdir in tranche_dirs:
            top_path = pipeline.rank_and_filter_tranche(tdir, filter_percent)
            self.assertEqual(pipeline.read_status(tdir), "DONE")
            self.assertTrue(top_path.exists())

            ranked_lines = (tdir / "results" / "ranked_all.tsv").read_text().splitlines()[1:]
            self.assertEqual(len(ranked_lines), 2)
            affinities = [float(line.split("\t")[1]) for line in ranked_lines]
            self.assertEqual(affinities, sorted(affinities))  # best (most negative) first

        combined_path = pipeline.combine_results(self.workdir, tranche_dirs)
        combined_lines = combined_path.read_text().splitlines()
        self.assertEqual(combined_lines[0], "ligand\taffinity_kcal_mol\ttranche")
        self.assertGreaterEqual(len(combined_lines), 2)

    def test_rerun_skips_completed_stages(self):
        setup, tranches = self._load()
        targets = io_parse.load_docking_targets(setup)
        tdir = pipeline.materialize_tranche(self.workdir, tranches[0], self.setup_path)

        pipeline.download_tranche(tdir)
        download_dir = tdir / "download"
        marker = next(download_dir.rglob("*.pdbqt.gz"))
        original_mtime = marker.stat().st_mtime

        # Re-running download on an already-DOWNLOADED tranche must be a no-op.
        pipeline.download_tranche(tdir)
        self.assertEqual(marker.stat().st_mtime, original_mtime)

        pipeline.unpack_tranche(tdir)
        pipeline.dock_tranche(
            tdir, targets, energy_range="3", cores=1, vinalc_bin=str(MOCK_VINALC), mpirun_bin=None
        )
        pipeline.rank_and_filter_tranche(tdir, float(setup["filter_percent"]))
        self.assertEqual(pipeline.read_status(tdir), "DONE")

        # Docking again on a DONE tranche should skip rather than re-invoke the binary.
        pipeline.dock_tranche(
            tdir, targets, energy_range="3", cores=1, vinalc_bin="/nonexistent/vinalc", mpirun_bin=None
        )
        self.assertEqual(pipeline.read_status(tdir), "DONE")

    def test_nuke_workdir_removes_directory(self):
        (self.workdir / "some_file.txt").write_text("hi")
        pipeline.nuke_workdir(self.workdir)
        self.assertFalse(self.workdir.exists())

    def test_nuke_workdir_missing_dir_is_noop(self):
        missing = self.tmp / "does_not_exist"
        pipeline.nuke_workdir(missing)  # should not raise

    def test_nuke_workdir_refuses_dangerous_paths(self):
        with self.assertRaises(RuntimeError):
            pipeline.nuke_workdir(Path.home())

    def test_clean_workdir_removes_only_tranches(self):
        setup, tranches = self._load()
        pipeline.materialize_tranche(self.workdir, tranches[0], self.setup_path)
        (self.workdir / "unrelated.txt").write_text("keep me")

        pipeline.clean_workdir(self.workdir)

        self.assertFalse((self.workdir / "tranches").exists())
        self.assertTrue((self.workdir / "unrelated.txt").exists())


if __name__ == "__main__":
    unittest.main()
