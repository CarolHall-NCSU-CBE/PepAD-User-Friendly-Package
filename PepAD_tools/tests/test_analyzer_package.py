import contextlib
import io
import inspect
import sys
import tempfile
import unittest
import warnings
from pathlib import Path
from unittest import mock

import matplotlib
import pandas as pd

matplotlib.use("Agg")

import pepad_tools.analyzer as analyzer_module

ROOT = Path(__file__).resolve().parents[1]


def load_python_functions() -> dict[str, object]:
    """Return functions from the installed analyzer package module."""
    return vars(analyzer_module)


def sample_energy_detail() -> pd.DataFrame:
    """Create accepted, MC-rejected, and rotamer-rejected trial records."""
    rows = [
        [1, 1, "trial", "Accept", "G-N-N", -2.0, -1.0, 0.2, 1.0, 0.5],
        [2, 1, "trial", "Accept", "G-N-Q", -1.5, 0.4, -0.2, 0.8, -0.3],
        [3, 1, "trial", "Reject-MC", "G-Q-N", 1.2, 1.0, -0.1, -0.5, -1.0],
        [4, 1, "trial", "Reject-MC", "Q-N-N", 1.8, -0.2, 0.3, -0.8, -0.6],
        [5, 1, "trial", "Reject-Rotamer", "Q-Q-N", None, None, None, None, None],
    ]
    columns = [
        "step",
        "attempt",
        "type",
        "acceptance",
        "Sequence",
        "dE_VDW",
        "dE_(ELE+SGB)",
        "dE_SUR",
        "dTS",
        "dPAGG",
    ]
    return pd.DataFrame(rows, columns=columns)


class PepADAnalyzerPythonTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.analyzer = load_python_functions()

    def test_public_functions_have_type_annotations(self) -> None:
        function_names = [
            "analyze_pdb",
            "read_input",
            "read_energy_profile",
            "plot_pepad",
            "generate_pepad_report",
            "get_parser",
            "read_energy_detail",
            "plot_energy_detail",
            "analyze_delta_energy",
            "main",
        ]

        for name in function_names:
            signature = inspect.signature(self.analyzer[name])
            self.assertIsNot(
                signature.return_annotation,
                inspect.Signature.empty,
                msg=f"{name} needs a return annotation",
            )
            for parameter in signature.parameters.values():
                self.assertIsNot(
                    parameter.annotation,
                    inspect.Signature.empty,
                    msg=f"{name}.{parameter.name} needs a type annotation",
                )

    def test_violin_argument_behavior(self) -> None:
        parser = self.analyzer["get_parser"]()

        self.assertIsNone(parser.parse_args(["--details"]).violin)
        self.assertEqual(parser.parse_args(["--details", "--violin"]).violin, [])
        self.assertEqual(
            parser.parse_args(
                ["--details", "--violin", "energy", "contribution"]
            ).violin,
            ["energy", "contribution"],
        )

    def test_bare_violin_selects_energy_plot(self) -> None:
        analyze = mock.Mock()
        originals = {
            "read_input": self.analyzer["read_input"],
            "read_energy_detail": self.analyzer["read_energy_detail"],
            "analyze_delta_energy": self.analyzer["analyze_delta_energy"],
        }
        self.analyzer["read_input"] = mock.Mock(
            return_value={
                "BASE_DIR": str(ROOT),
                "N_AA": 3,
                "N_RESIDUES": 3,
            }
        )
        self.analyzer["read_energy_detail"] = mock.Mock(
            return_value=sample_energy_detail()
        )
        self.analyzer["analyze_delta_energy"] = analyze

        try:
            with mock.patch.object(
                sys, "argv", ["analyzer", "--details", "--violin"]
            ):
                self.analyzer["main"]()
        finally:
            self.analyzer.update(originals)

        self.assertEqual(analyze.call_args.args[2], ["energy"])

    def test_top_argument_reaches_profile_report(self) -> None:
        profile = pd.DataFrame({"Sequence": ["GNN"], "Score": [-1.0]})
        parameters = {
            "BASE_DIR": str(ROOT),
            "N_AA": 3,
            "N_RESIDUES": 3,
        }
        report = mock.Mock()
        originals = {
            "read_input": self.analyzer["read_input"],
            "read_energy_profile": self.analyzer["read_energy_profile"],
            "generate_pepad_report": self.analyzer["generate_pepad_report"],
        }
        self.analyzer["read_input"] = mock.Mock(return_value=parameters)
        self.analyzer["read_energy_profile"] = mock.Mock(return_value=profile)
        self.analyzer["generate_pepad_report"] = report

        try:
            with mock.patch.object(
                sys, "argv", ["analyzer", "--profiles", "--top", "200"]
            ):
                self.analyzer["main"]()
        finally:
            self.analyzer.update(originals)

        self.assertIs(report.call_args.args[0], profile)
        self.assertIs(report.call_args.args[1], parameters)
        self.assertEqual(report.call_args.args[2], 200)

    def test_profile_report_uses_underscore_filename(self) -> None:
        profile = pd.DataFrame(
            {
                "Sequence": ["GNN", "GNN", "NNG"],
                "Score": [-1.0, -2.0, -0.5],
            }
        )
        with tempfile.TemporaryDirectory() as temp_dir:
            parameters = {"BASE_DIR": temp_dir}
            self.analyzer["generate_pepad_report"](
                profile, parameters, top=2
            )

            self.assertTrue(
                (Path(temp_dir) / "PepAD_report.txt").is_file()
            )
            self.assertFalse(
                (Path(temp_dir) / "PepAD report.txt").exists()
            )

    def test_read_energy_detail_includes_rotamer_rejection(self) -> None:
        content = (
            "1 1 -1 1 2 3 4 5 6 -1 -2 -3 -4 5 6 trial Accept\n"
            "G N N Q Q N Y\n"
            "2 1 trial Reject-Rotamer\n"
            "G N N Q Q N Y\n"
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            detail_file = Path(temp_dir) / "energydetails.txt"
            detail_file.write_text(content, encoding="utf-8")
            result = self.analyzer["read_energy_detail"](temp_dir)

        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(
            result["acceptance"].tolist(), ["Accept", "Reject-Rotamer"]
        )
        self.assertTrue(pd.isna(result.loc[1, "dE_VDW"]))

    def test_energy_profile_removes_incomplete_records(self) -> None:
        content = (
            "0 GNN 0 0 0 0\n"
            "1 GNN ****** ****** -5 0.0\n"
            "2 GNN -1 -2 -3 0.1\n"
            "3 GNN -1 -2\n"
            "4 GNN -1 -2 -3 0.2\n"
            "5 GNU -1 -2 -3 0.3\n"
            "****** GNN -1 -2 -3 0.4\n"
            "6 GNN -1 -2 -3 0.5\n"
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            profile_file = Path(temp_dir) / "energyprofile.txt"
            profile_file.write_text(content, encoding="utf-8")
            output = io.StringIO()
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                with contextlib.redirect_stdout(output):
                    profile = self.analyzer["read_energy_profile"](temp_dir, 3)

        self.assertEqual(profile["step"].tolist(), [0, 1, 2, 3])
        self.assertEqual(len(caught), 1)
        warning = str(caught[0].message)
        self.assertIn("file line 2, step 1", warning)
        self.assertIn("invalid numeric value(s) in Score, G_bind", warning)
        self.assertIn("file line 4, step 3", warning)
        self.assertIn("file line 6, step 5", warning)
        self.assertIn("file line 7, step unknown", warning)
        self.assertIn("Removed records are treated as missing steps", warning)

    def test_input_warns_when_n_aa_disagrees_with_pdb(self) -> None:
        input_text = (
            "PDBFILE = peptide.pdb\n"
            "N_STEPS = 100\n"
            "N_AA = 8\n"
        )
        original_analyze_pdb = self.analyzer["analyze_pdb"]
        self.analyzer["analyze_pdb"] = mock.Mock(return_value=(7, 9, 4))

        try:
            with tempfile.TemporaryDirectory() as temp_dir:
                input_file = Path(temp_dir) / "input.txt"
                input_file.write_text(input_text, encoding="utf-8")
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always")
                    parameters = self.analyzer["read_input"](temp_dir)
        finally:
            self.analyzer["analyze_pdb"] = original_analyze_pdb

        self.assertEqual(len(caught), 1)
        self.assertIn("does not match", str(caught[0].message))
        self.assertEqual(parameters["N_AA"], 7)
        self.assertEqual(parameters["N_RESIDUES"], 9)
        self.assertEqual(parameters["N_CHAINS"], 4)

    def test_energy_detail_rejects_incomplete_records(self) -> None:
        full_energy = (
            "1 1 -1 1 2 3 4 5 6 -1 -2 -3 -4 5 6 trial Accept\n"
        )
        cases = {
            "missing energy terms": (
                "2 4 -1.1519 -17.6428 -43.9508 40.7703\n"
                "G N N Q Q N Y\n",
                "Missing or extra value at line 1",
            ),
            "missing sequence": (
                full_energy,
                "Each record must contain an energy line and a sequence line",
            ),
            "partial sequence": (
                full_energy + "G N N Q Q N\n",
                "has 6 residues; expected 7",
            ),
        }

        for name, (content, message) in cases.items():
            with self.subTest(name=name):
                with tempfile.TemporaryDirectory() as temp_dir:
                    detail_file = Path(temp_dir) / "energydetails.txt"
                    detail_file.write_text(content, encoding="utf-8")
                    with self.assertRaisesRegex(ValueError, message):
                        self.analyzer["read_energy_detail"](temp_dir, 7)

    def test_profile_steps_are_normalized(self) -> None:
        content = (
            "0 GNN 0 0 0 0\n"
            "1 GNN 0 0 0 0\n"
            "1 GNQ 0 0 0 0\n"
            "2 GNQ 0 0 0 0\n"
            "5 GQN 0 0 0 0\n"
            "6 GQN 0 0 0 0\n"
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            profile_file = Path(temp_dir) / "energyprofile.txt"
            profile_file.write_text(content, encoding="utf-8")
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                profile = self.analyzer["read_energy_profile"](temp_dir, 3)

        self.assertEqual(profile["step"].tolist(), [0, 1, 2, 3, 4, 5])
        self.assertIn("Duplicate step 1 found", output.getvalue())
        self.assertIn("shifted back by 1", output.getvalue())
        self.assertIn("Missing step(s) 3-4 found", output.getvalue())
        self.assertIn("shifted forward by 2", output.getvalue())
        self.assertIn("Corrected final step: 5", output.getvalue())

    def test_energy_detail_step_blocks_are_normalized(self) -> None:
        def record(step: int, attempt: int) -> str:
            energy = (
                f"{step} {attempt} -1 1 2 3 4 5 6 "
                "-1 -2 -3 -4 5 6 trial Accept\n"
            )
            return energy + "G N N\n"

        content = "".join(
            [
                record(1, 1),
                record(1, 2),
                record(2, 1),
                record(2, 2),
                record(2, 1),
                record(2, 2),
                record(3, 1),
                record(6, 1),
                record(7, 1),
            ]
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            detail_file = Path(temp_dir) / "energydetails.txt"
            detail_file.write_text(content, encoding="utf-8")
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                detail = self.analyzer["read_energy_detail"](temp_dir, 3)

        self.assertEqual(
            detail["step"].tolist(), [1, 1, 2, 2, 3, 3, 4, 5, 6]
        )
        self.assertIn("Duplicate step 2 found", output.getvalue())
        self.assertIn("Missing step(s) 4-5 found", output.getvalue())
        self.assertIn("Corrected final step: 6", output.getvalue())

    def test_details_always_write_report_without_plots(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            report = Path(temp_dir) / "Detail_report.txt"
            result = self.analyzer["analyze_delta_energy"](
                sample_energy_detail(), str(report), []
            )

            self.assertTrue(report.is_file())
            self.assertFalse(
                (Path(temp_dir) / "Delta_energy_distribution.png").exists()
            )
            self.assertFalse(
                (Path(temp_dir) / "Delta_energy_contribution.png").exists()
            )
            expected_header = (
                f"{'Ei':<9}{'Median_ddEi(-)(kcal/mol)':>25}"
                f"{'Mdn_ddEi(+)(kcal/mol)':>25}"
                f"{'Frequency_of_ddEi(-)(%)':>25}"
                f"{'Mdn_[|ddEi(-)|/sum_i|ddEi(-)|](%)':>35}"
                f"{'Frequency_of_ddEi(+)(%)':>25}"
                f"{'Mdn_[|ddEi(+)|/sum_i|ddEi(+)|](%)':>35}"
            )
            self.assertIn(
                expected_header, report.read_text(encoding="utf-8")
            )
            self.assertNotIn(
                "Reject-MC trials without a positive contribution",
                report.read_text(encoding="utf-8"),
            )

        self.assertEqual(len(result), 3)
        self.assertTrue(all(isinstance(item, pd.DataFrame) for item in result))

    def test_violin_order_and_600_dpi(self) -> None:
        figure_type = matplotlib.figure.Figure

        with tempfile.TemporaryDirectory() as temp_dir:
            report = Path(temp_dir) / "Detail_report.txt"
            with mock.patch.object(figure_type, "savefig") as savefig:
                with mock.patch.object(self.analyzer["plt"], "show"):
                    self.analyzer["analyze_delta_energy"](
                        sample_energy_detail(),
                        str(report),
                        ["energy", "contribution"],
                    )

        filenames = [
            Path(call.args[0]).name for call in savefig.call_args_list
        ]
        self.assertEqual(
            filenames,
            [
                "Delta_energy_distribution.png",
                "Delta_energy_contribution.png",
            ],
        )
        self.assertTrue(
            all(call.kwargs["dpi"] == 600 for call in savefig.call_args_list)
        )

    def test_missing_energy_detail_raises_clear_error(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            with self.assertRaisesRegex(
                FileNotFoundError, "energydetails.txt is not found"
            ):
                self.analyzer["read_energy_detail"](temp_dir)


if __name__ == "__main__":
    unittest.main()
