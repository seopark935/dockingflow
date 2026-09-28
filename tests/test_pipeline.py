#!/usr/bin/env python3
"""End-to-end tests for dockingflow, run entirely offline with a mock vinalc.

Run with: python3 -m unittest discover -s tests -v
(from the repo root)
"""
from __future__ import annotations

import gzip
import json
import shutil
import threading
import urllib.error
import urllib.request
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

import io_parse
import pipeline
import resources
import web_server
import zinc_split

MOCK_VINALC = REPO_ROOT / "tests" / "fixtures" / "mock_vinalc.py"
# A real ZINC22 3D archive (4 molecules), downloaded from
# https://files.docking.org/zinc22/zinc-22a/H04/H04M000/a/H04M000-N-aaaaaa.pdbqt.tgz
ZINC22_SAMPLE = REPO_ROOT / "tests" / "fixtures" / "H04M000-N-aaaaaa.pdbqt.tgz"


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

    def _dock(self, tdir, targets, setup, vinalc_bin=str(MOCK_VINALC)):
        pipeline.dock_tranche(
            tdir, targets, options=io_parse.load_vinalc_options(setup), vinalc_bin=vinalc_bin, mpirun_bin=None
        )

    def test_validate_inputs_accepts_repo_fixtures(self):
        setup, tranches = self._load()
        self.assertEqual(len(tranches), 1)
        self.assertEqual(tranches[0].log_p, 5.0)
        self.assertEqual(tranches[0].molecular_weight, 200)

    def test_load_setup_resolves_paths_relative_to_setup_file(self):
        setup = io_parse.load_setup(self.setup_path)
        self.assertEqual(Path(setup["recList"]), REPO_ROOT / "recList.txt")

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
            self.assertGreater(pipeline.count_ligand_archives(tdir / "download"), 0)

        for tdir in tranche_dirs:
            pipeline.unpack_tranche(tdir)
            self.assertEqual(pipeline.read_status(tdir), "UNPACKED")
            self.assertEqual(len((tdir / "ligand_list.txt").read_text().splitlines()), 1)
            index = pipeline.read_ligand_index(tdir)
            self.assertEqual(index, {1: "ZINC450000002E48.0", 2: "ZINC450000002Jfn.0"})

        for tdir in tranche_dirs:
            self._dock(tdir, targets, setup)
            self.assertEqual(pipeline.read_status(tdir), "DOCKED")
            self.assertTrue((tdir / "docking" / pipeline.VINALC_POSES).exists())

        filter_percent = float(setup["filter_percent"])
        for tdir in tranche_dirs:
            top_path = pipeline.rank_and_filter_tranche(tdir, filter_percent)
            self.assertEqual(pipeline.read_status(tdir), "DONE")
            self.assertTrue(top_path.exists())

            ranked_lines = (tdir / "results" / "ranked_all.tsv").read_text().splitlines()[1:]
            self.assertEqual({line.split("\t")[0] for line in ranked_lines}, {"ZINC450000002E48.0", "ZINC450000002Jfn.0"})
            affinities = [float(line.split("\t")[1]) for line in ranked_lines]
            self.assertEqual(affinities, sorted(affinities))  # best (most negative) first
            self.assertTrue(all(line.split("\t")[2] == "protein" for line in ranked_lines))

        combined_path = pipeline.combine_results(self.workdir, tranche_dirs)
        combined_lines = combined_path.read_text().splitlines()
        self.assertEqual(combined_lines[0], "ligand\taffinity_kcal_mol\treceptor\ttranche")
        self.assertGreaterEqual(len(combined_lines), 2)

    def test_rerun_skips_completed_stages(self):
        setup, tranches = self._load()
        targets = io_parse.load_docking_targets(setup)
        tdir = pipeline.materialize_tranche(self.workdir, tranches[0], self.setup_path)

        pipeline.download_tranche(tdir)
        marker = pipeline.find_ligand_archives(tdir / "download")[0]
        original_mtime = marker.stat().st_mtime

        # Re-running download on an already-DOWNLOADED tranche must be a no-op.
        pipeline.download_tranche(tdir)
        self.assertEqual(marker.stat().st_mtime, original_mtime)

        pipeline.unpack_tranche(tdir)
        self._dock(tdir, targets, setup)
        pipeline.rank_and_filter_tranche(tdir, float(setup["filter_percent"]))
        self.assertEqual(pipeline.read_status(tdir), "DONE")

        # Docking again on a DONE tranche should skip rather than re-invoke the binary.
        self._dock(tdir, targets, setup, vinalc_bin="/nonexistent/vinalc")
        self.assertEqual(pipeline.read_status(tdir), "DONE")

    def test_partial_download_failure_is_recorded_not_fatal(self):
        tdir = self.tmp / "tranche"
        (tdir / "inputs").mkdir(parents=True)
        (tdir / "inputs" / "curl_script.curl").write_text(
            zinc_split.SCRIPT_HEADER.format(source="test", code="H04M000")
            + "echo hi > ok.pdbqt\nfalse\n"
            + zinc_split.SCRIPT_FOOTER
        )
        pipeline.download_tranche(tdir)
        self.assertEqual(pipeline.read_status(tdir), "DOWNLOADED")
        self.assertIn("FAILED", (tdir / "download_failures.txt").read_text())

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


class UnpackFormatTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="dockingflow_unpack_"))
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def _unpack(self, files: dict[str, Path | bytes]) -> Path:
        tdir = self.tmp / "tranche"
        download = tdir / "download"
        download.mkdir(parents=True)
        for name, src in files.items():
            dest = download / name
            dest.parent.mkdir(parents=True, exist_ok=True)
            if isinstance(src, Path):
                shutil.copy(src, dest)
            else:
                dest.write_bytes(src)
        pipeline.write_status(tdir, "DOWNLOADED")
        pipeline.unpack_tranche(tdir)
        return tdir

    def test_real_zinc22_archive_is_wrapped_in_models(self):
        tdir = self._unpack({"H04/H04M000/a/H04M000-N-aaaaaa.pdbqt.tgz": ZINC22_SAMPLE})
        index = pipeline.read_ligand_index(tdir)
        self.assertEqual(len(index), 4)
        self.assertTrue(all(name.startswith("ZINC45") for name in index.values()))

        # VinaLC only docks MODEL...ENDMDL blocks: every molecule must be wrapped.
        text = (tdir / "ligands" / "ligands_00001.pdbqt").read_text()
        self.assertEqual(text.count("\nENDMDL"), 4)
        self.assertEqual(sum(1 for ln in text.splitlines() if ln.startswith("MODEL")), 4)

    def test_multi_model_pdbqt_gz_and_corrupt_archive(self):
        multi = (
            "MODEL 1\nREMARK  Name = ZINC000000000001\nROOT\nATOM      1  C   LIG    1       0.0 0.0 0.0  0.00  0.00    +0.000 C\nENDROOT\nTORSDOF 0\nENDMDL\n"
            "MODEL 2\nROOT\nATOM      1  C   LIG    1       1.0 1.0 1.0  0.00  0.00    +0.000 C\nENDROOT\nTORSDOF 0\nENDMDL\n"
        )
        tdir = self._unpack({
            "AA/AAAA.xaa.pdbqt.gz": gzip.compress(multi.encode()),
            "broken.pdbqt.tgz": b"<html>404 Not Found</html>",
        })
        index = pipeline.read_ligand_index(tdir)
        self.assertEqual(index, {1: "ZINC000000000001", 2: "AAAA.xaa_2"})
        self.assertIn("broken.pdbqt.tgz", (tdir / "unpack_warnings.txt").read_text())
        # An unnamed molecule gets a Name remark so VinaLC's output poses carry its id.
        self.assertIn("REMARK  Name = AAAA.xaa_2", (tdir / "ligands" / "ligands_00001.pdbqt").read_text())


class VinaLCOutputTest(unittest.TestCase):
    def test_parse_poses_maps_out_of_order_records(self):
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp, ignore_errors=True)
        poses = tmp / "out.pdbqt.gz"
        with gzip.open(poses, "wt") as f:
            f.write(
                "REMARK RECEPTOR /x/recB.pdbqt\nREMARK LIGAND LIGAND 2\n"
                "MODEL 1\nREMARK VINA RESULT:      -9.1      0.000      0.000\nENDMDL\n"
                "MODEL 2\nREMARK VINA RESULT:      -8.0      1.000      2.000\nENDMDL\n\n"
                "REMARK RECEPTOR /x/recA.pdbqt\nREMARK LIGAND LIGAND 1\n\n"  # no pose found
                "REMARK RECEPTOR /x/recA.pdbqt\nREMARK LIGAND LIGAND 2\n"
                "MODEL 1\nREMARK VINA RESULT:      -6.5      0.000      0.000\nENDMDL\n"
            )
        self.assertEqual(
            list(pipeline.parse_vinalc_poses(poses)),
            [("/x/recB.pdbqt", 2, -9.1), ("/x/recA.pdbqt", 2, -6.5)],
        )

    def test_command_uses_real_vinalc_flags_and_rank_count(self):
        opts = io_parse.VinaLCOptions(mpi_ranks=25, energy_range="3", exhaustiveness=8, num_modes=9, seed="42")
        cmd = pipeline.build_vinalc_command(opts, "vinalc", "mpirun")
        self.assertEqual(cmd[:3], ["mpirun", "-np", "25"])
        for flag in ("--recList", "--ligList", "--geoList", "--exhaustiveness", "--energy_range", "--seed"):
            self.assertIn(flag, cmd)
        self.assertNotIn("--config", cmd)


class ResourcePlanningTest(unittest.TestCase):
    def test_ranks_sized_by_exhaustiveness_not_one_per_core(self):
        self.assertEqual(io_parse.plan_mpi_ranks(192, 8, None, 1.0), 25)
        self.assertEqual(io_parse.plan_mpi_ranks(4, 8, None, 1.0), 2)  # VinaLC needs >= 2 ranks

    def test_memory_budget_caps_ranks(self):
        self.assertEqual(io_parse.plan_mpi_ranks(192, 8, 10.2, 1.0), 11)

    def test_setup_memory_gb_is_applied(self):
        setup = io_parse.load_setup(REPO_ROOT / "setup.txt")
        uncapped = io_parse.load_vinalc_options(setup).mpi_ranks
        setup["memory_gb"] = "4"
        self.assertLess(io_parse.load_vinalc_options(setup).mpi_ranks, uncapped)

    def test_parse_server_report_and_recommend(self):
        report = (
            "===== DOCKINGFLOW SUMMARY =====\nDF_HOSTNAME=dock1\nDF_LOGICAL_CPUS=384\nDF_PHYSICAL_CORES=192\n"
            "DF_MEM_TOTAL_KB=1048576000\nDF_MEM_AVAILABLE_KB=943718400\nDF_LOAD1=4.0\n"
        )
        res = resources.parse_server_report(report)
        self.assertEqual((res["hostname"], res["physical_cores"]), ("dock1", 192))
        self.assertAlmostEqual(res["mem_total_gb"], 1000.0)
        rec = resources.recommend(res)
        self.assertEqual(rec["cores"], 192 - 4 - 9)
        self.assertAlmostEqual(rec["memory_gb"], 720.0)
        with self.assertRaises(ValueError):
            resources.parse_server_report("not a report")

    def test_write_budget_to_setup_keeps_other_lines(self):
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp, ignore_errors=True)
        setup = tmp / "setup.txt"
        setup.write_text("# comment\nrecList=recList.txt\ncores=192\n")
        resources.write_budget_to_setup(str(setup), 64, 128.0)
        self.assertEqual(setup.read_text(), "# comment\nrecList=recList.txt\ncores=64\nmemory_gb=128\n")


class ZincSplitTest(unittest.TestCase):
    def test_groups_commands_by_tranche_and_writes_loadable_map(self):
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp, ignore_errors=True)
        url = "https://files.docking.org/zinc22/zinc-22a"
        downloader = tmp / "ZINC22-downloader-3D-pdbqt.tgz.curl"
        downloader.write_text(
            f"curl --fail --create-dirs -o H04/H04M000/a/H04M000-N-aaaaaa.pdbqt.tgz {url}/H04/H04M000/a/H04M000-N-aaaaaa.pdbqt.tgz\n"
            f"curl --fail --create-dirs -o H04/H04M000/a/H04M000-O-aaaaaa.pdbqt.tgz {url}/H04/H04M000/a/H04M000-O-aaaaaa.pdbqt.tgz\n"
            f"curl --fail --create-dirs -o H17/H17P050/a/H17P050-N-aaaaaa.pdbqt.tgz {url}/H17/H17P050/a/H17P050-N-aaaaaa.pdbqt.tgz\n"
        )
        subprocess.run(
            [sys.executable, str(REPO_ROOT / "zinc_split.py"), str(downloader),
             "--out-dir", str(tmp / "scripts"), "--map", str(tmp / "tranches.txt")],
            check=True, capture_output=True,
        )
        tranches = io_parse.load_tranches_tsv(tmp / "tranches.txt")
        self.assertEqual([io_parse.tranche_label(t) for t in tranches], ["H04M000", "H17P050"])
        self.assertEqual((tranches[1].heavy_atoms, tranches[1].log_p), (17, 0.5))
        self.assertEqual((tmp / "scripts" / "H04M000.curl").read_text().count("curl --fail"), 2)


class WebServerTest(unittest.TestCase):
    """The HTTP layer behind `gui.py --web` (the GUI on a headless server)."""

    class FakeAPI:
        def get_defaults(self):
            return {"workdir": "/x"}

        def add(self, a, b):
            return a + b

        def nuke(self, workdir):
            return {"ok": True}

        def set_window(self, w):
            return None

        def _run(self):
            return None

    def setUp(self):
        from http.server import ThreadingHTTPServer

        self.token = "secret-token"
        handler = web_server.make_handler(self.FakeAPI(), REPO_ROOT / "gui_assets" / "index.html", self.token)
        self.server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
        threading.Thread(target=self.server.serve_forever, daemon=True).start()
        self.addCleanup(self.server.server_close)
        self.addCleanup(self.server.shutdown)
        self.base = f"http://127.0.0.1:{self.server.server_address[1]}"

    def _post(self, method, args, token="secret-token"):
        req = urllib.request.Request(
            f"{self.base}/api/{method}", data=json.dumps(args).encode(), method="POST",
            headers={web_server.TOKEN_HEADER: token},
        )
        try:
            with urllib.request.urlopen(req) as r:
                return r.status, json.load(r)
        except urllib.error.HTTPError as e:
            with e:
                return e.code, json.load(e)

    def test_page_is_served_in_web_mode(self):
        with urllib.request.urlopen(self.base + "/?token=abc") as r:
            self.assertIn("window.DF_WEB = true;</script>", r.read().decode())

    def test_api_call_with_token(self):
        self.assertEqual(self._post("add", [2, 3]), (200, 5))
        self.assertEqual(self._post("get_defaults", []), (200, {"workdir": "/x"}))

    def test_api_rejects_missing_or_wrong_token(self):
        self.assertEqual(self._post("nuke", ["/x"], token="")[0], 403)
        self.assertEqual(self._post("nuke", ["/x"], token="guess")[0], 403)

    def test_private_and_desktop_only_methods_not_exposed(self):
        self.assertEqual(self._post("_run", [])[0], 404)
        self.assertEqual(self._post("set_window", [None])[0], 404)

    def test_list_dir_for_in_page_picker(self):
        import gui

        res = gui.PipelineAPI().list_dir(str(REPO_ROOT / "setup.txt"))  # a file lists its directory
        self.assertTrue(res["ok"])
        self.assertEqual(Path(res["path"]), REPO_ROOT)
        names = [e["name"] for e in res["entries"]]
        self.assertIn("setup.txt", names)
        self.assertTrue(res["entries"][0]["is_dir"])  # directories first


if __name__ == "__main__":
    unittest.main()
