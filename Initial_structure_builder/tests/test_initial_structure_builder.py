import contextlib
import functools
import importlib.util
import io
import inspect
import os
import random
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
BUILDER_PATH = PROJECT_ROOT / "initial_structure_builder.py"
ROTAMER_LIBRARY = PROJECT_ROOT / "RotamerLibrary"
ATOM_COLUMNS = {"anum", "atom_name", "aa_name", "resid", "x", "y", "z"}


def load_builder_functions():
    """Load initial_structure_builder.py without running its terminal interface."""
    if not BUILDER_PATH.is_file():
        raise FileNotFoundError(
            f"Initial structure builder does not exist: {BUILDER_PATH}"
        )

    module_name = "initial_structure_builder_under_test"
    spec = importlib.util.spec_from_file_location(module_name, BUILDER_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load builder Python file: {BUILDER_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    namespace = vars(module)

    original_builder = namespace[
        "build_single_peptide_with_local_peptidebuilder"
    ]

    @functools.wraps(original_builder)
    def build_with_project_rotamers(sequence, angles):
        return original_builder(
            sequence,
            angles,
            rotamer_dir=ROTAMER_LIBRARY,
        )

    namespace["build_single_peptide_with_local_peptidebuilder"] = (
        build_with_project_rotamers
    )
    return namespace


class InitialStructureBuilderTests(unittest.TestCase):
    """Check the main builder inputs, outputs, and generated files."""

    @classmethod
    def setUpClass(cls):
        cls.builder = load_builder_functions()

    def setUp(self):
        random.seed(1234)
        self.original_directory = Path.cwd()
        self.temp_directory = tempfile.TemporaryDirectory(
            prefix="initial_structure_builder_"
        )
        os.chdir(self.temp_directory.name)

    def tearDown(self):
        os.chdir(self.original_directory)
        self.temp_directory.cleanup()

    def assert_atom_dataframe(self, atoms):
        pd = self.builder["pd"]
        np = self.builder["np"]
        self.assertIsInstance(atoms, pd.DataFrame)
        self.assertFalse(atoms.empty)
        self.assertTrue(ATOM_COLUMNS.issubset(atoms.columns))
        coordinates = atoms[["x", "y", "z"]].to_numpy(dtype=float)
        self.assertTrue(np.isfinite(coordinates).all())

    def assert_terminal_style(self, atoms, cap_flag):
        residue_names = set(atoms["aa_name"].astype(str))
        if cap_flag == 0:
            self.assertTrue({"ACE", "NME", "NHE"}.isdisjoint(residue_names))
        elif cap_flag == 1:
            self.assertTrue({"ACE", "NME"}.issubset(residue_names))
            self.assertNotIn("NHE", residue_names)
        else:
            self.assertTrue({"ACE", "NHE"}.issubset(residue_names))
            self.assertNotIn("NME", residue_names)

    def run_build_case(self, case_name, output_name, **builder_inputs):
        case_directory = Path(self.temp_directory.name) / case_name
        case_directory.mkdir()
        previous_directory = Path.cwd()
        output = io.StringIO()
        try:
            os.chdir(case_directory)
            with contextlib.redirect_stdout(output):
                atoms = self.builder["build_sheets"](
                    "GNNQQNY",
                    output_name,
                    chains=2,
                    **builder_inputs,
                )
        finally:
            os.chdir(previous_directory)

        pdb_files = {path.name for path in case_directory.glob("*.pdb")}
        return atoms, pdb_files, output.getvalue()

    def test_class1_build_returns_atoms_and_expected_files(self):
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            atoms = self.builder["build_sheets"](
                "GNNQQNY",
                "class1",
                classes=1,
                chains=2,
                cap_flag=0,
            )

        self.assert_atom_dataframe(atoms)
        self.assertEqual(
            {path.name for path in Path.cwd().glob("*.pdb")},
            {"para_A.pdb", "class1.pdb"},
        )
        summary = output.getvalue()
        self.assertIn("Sheets structure: Class 1", summary)
        self.assertIn('Sequence: "GNNQQNY"', summary)
        self.assertIn("Strands per sheet: 2", summary)

    def test_two_sequence_hybrid_writes_only_required_pdb_files(self):
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            atoms = self.builder["build_sheets"](
                "GNNQQNY",
                "hybrid1",
                sequence_b="AAAAAAA",
                pattern1="AB",
                pattern2="BA",
                hybrid_type=1,
                cap_flag=1,
            )

        self.assert_atom_dataframe(atoms)
        self.assertEqual(
            {path.name for path in Path.cwd().glob("*.pdb")},
            {
                "para_A.pdb",
                "anti_A.pdb",
                "para_B.pdb",
                "anti_B.pdb",
                "hybrid1.pdb",
            },
        )
        summary = output.getvalue()
        self.assertIn("Sheets structure: Hybrid 1", summary)
        self.assertIn('Sequence A: "GNNQQNY"', summary)
        self.assertIn('Sequence B: "AAAAAAA"', summary)
        self.assertNotIn("temp_", " ".join(path.name for path in Path.cwd().iterdir()))

    def test_terminal_processing_stays_in_memory(self):
        build_single = self.builder[
            "build_single_peptide_with_local_peptidebuilder"
        ]
        peptide = build_single("GNNQQNY", [-119, 113])
        uncapped = self.builder["add_NH_remove_OH"](peptide)
        ace_nme = self.builder["add_caps"](peptide, 1)
        ace_nhe = self.builder["add_caps"](peptide, 2)

        for atoms in (peptide, uncapped, ace_nme, ace_nhe):
            self.assert_atom_dataframe(atoms)
        self.assertTrue(uncapped["aa_name"].astype(str).str.startswith("N").any())
        self.assertTrue(uncapped["aa_name"].astype(str).str.startswith("C").any())
        self.assertTrue({"ACE", "NME"}.issubset(set(ace_nme["aa_name"])))
        self.assertTrue({"ACE", "NHE"}.issubset(set(ace_nhe["aa_name"])))
        self.assertEqual(list(Path.cwd().glob("*.pdb")), [])

    def test_invalid_inputs_raise_clear_errors(self):
        build_sheets = self.builder["build_sheets"]
        with self.assertRaisesRegex(ValueError, "invalid one-letter"):
            build_sheets("GX", "invalid", classes=1)
        with self.assertRaisesRegex(ValueError, "requires both pattern1 and pattern2"):
            build_sheets(
                "GNNQQNY",
                "invalid",
                sequence_b="AAAAAAA",
                classes=1,
            )
        with self.assertRaisesRegex(ValueError, "either classes or hybrid_type"):
            build_sheets(
                "GNNQQNY",
                "invalid",
                classes=1,
                hybrid_type=1,
            )

    def test_missing_rotamer_paths_raise_clear_errors(self):
        read_rotamer = self.builder["read_rotamer_template"]

        missing_library = Path(self.temp_directory.name) / "missing_library"
        with self.assertRaises(FileNotFoundError) as folder_error:
            read_rotamer("ALA", missing_library)
        folder_message = str(folder_error.exception)
        self.assertIn("RotamerLibrary folder does not exist", folder_message)
        self.assertIn(str(missing_library.resolve()), folder_message)
        self.assertIn("entire RotamerLibrary folder", folder_message)

        empty_library = Path(self.temp_directory.name) / "empty_library"
        empty_library.mkdir()
        missing_file = empty_library / "ALA"
        with self.assertRaises(FileNotFoundError) as file_error:
            read_rotamer("ALA", empty_library)
        file_message = str(file_error.exception)
        self.assertIn("Rotamer file does not exist", file_message)
        self.assertIn(str(missing_file.resolve()), file_message)
        self.assertIn("residue files inside RotamerLibrary", file_message)

    def test_core_functions_have_input_and_return_annotations(self):
        function_names = (
            "add_NH_remove_OH",
            "add_caps",
            "read_rotamer_template",
            "build_target_beta_backbone",
            "Add_backbone_hydrogen",
            "Transplant",
            "build_single_peptide_with_local_peptidebuilder",
            "alignment",
            "_build_aligned_peptide_template",
            "_finite_float",
            "validate_sheet_inputs",
            "build_sheets",
            "read_arguments",
        )

        for name in function_names:
            with self.subTest(function=name):
                signature = inspect.signature(self.builder[name])
                for parameter in signature.parameters.values():
                    self.assertIsNot(
                        parameter.annotation,
                        inspect.Parameter.empty,
                        f"{name}.{parameter.name} has no type annotation",
                    )
                self.assertIsNot(
                    signature.return_annotation,
                    inspect.Signature.empty,
                    f"{name} has no return annotation",
                )

    def test_confirmed_legacy_functions_are_removed(self):
        removed_functions = (
            "linear_fitting_3D_points",
            "_write_simple_pdb",
            "copy_pdb_file",
            "packpep",
            "packpep_antiparallel",
        )
        for name in removed_functions:
            self.assertNotIn(name, self.builder)

    def test_all_classes_hybrids_and_terminal_styles(self):
        parallel_classes = {1, 2, 3, 4, 9, 10}

        for cap_flag in (0, 1, 2):
            for classes in range(1, 11):
                case_name = f"class_{classes}_cap_{cap_flag}"
                with self.subTest(classes=classes, cap_flag=cap_flag):
                    random.seed(case_name)
                    atoms, pdb_files, summary = self.run_build_case(
                        case_name,
                        case_name,
                        classes=classes,
                        cap_flag=cap_flag,
                    )
                    template = (
                        "para_A.pdb"
                        if classes in parallel_classes
                        else "anti_A.pdb"
                    )
                    self.assert_atom_dataframe(atoms)
                    self.assert_terminal_style(atoms, cap_flag)
                    self.assertEqual(
                        pdb_files,
                        {template, f"{case_name}.pdb"},
                    )
                    self.assertIn(
                        f"Sheets structure: Class {classes}",
                        summary,
                    )

            for hybrid_type in range(1, 7):
                case_name = f"hybrid_{hybrid_type}_cap_{cap_flag}"
                with self.subTest(
                    hybrid_type=hybrid_type,
                    cap_flag=cap_flag,
                ):
                    random.seed(case_name)
                    atoms, pdb_files, summary = self.run_build_case(
                        case_name,
                        case_name,
                        hybrid_type=hybrid_type,
                        cap_flag=cap_flag,
                    )
                    self.assert_atom_dataframe(atoms)
                    self.assert_terminal_style(atoms, cap_flag)
                    self.assertEqual(
                        pdb_files,
                        {
                            "para_A.pdb",
                            "anti_A.pdb",
                            f"{case_name}.pdb",
                        },
                    )
                    self.assertIn(
                        f"Sheets structure: Hybrid {hybrid_type}",
                        summary,
                    )


if __name__ == "__main__":
    unittest.main()
