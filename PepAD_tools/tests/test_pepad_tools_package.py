import contextlib
import io
import os
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

import pepad_tools
import pepad_tools.builder as builder_module


class PepadToolsPackageTests(unittest.TestCase):
    """Check installed-package resources, imports, and structure building."""

    def setUp(self):
        self.original_directory = Path.cwd()
        self.temp_directory = tempfile.TemporaryDirectory(
            prefix="pepad_tools_package_"
        )
        os.chdir(self.temp_directory.name)

    def tearDown(self):
        os.chdir(self.original_directory)
        self.temp_directory.cleanup()

    def test_package_metadata_and_rotamer_data(self):
        self.assertEqual(pepad_tools.__version__, "1.42")
        self.assertTrue(builder_module.DEFAULT_ROTAMER_DIR.is_dir())
        self.assertEqual(
            len(list(builder_module.DEFAULT_ROTAMER_DIR.iterdir())),
            30,
        )

        alanine = builder_module.read_rotamer_template("ALA")
        self.assertIsInstance(alanine, pd.DataFrame)
        self.assertFalse(alanine.empty)

    def test_help_names_pepad_and_amber_formats(self):
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            with self.assertRaises(SystemExit) as exit_context:
                builder_module.read_arguments(["--help"])

        self.assertEqual(exit_context.exception.code, 0)
        self.assertIn(
            "PDB format: 0 = PepAD format; 1 = AMBER format.",
            output.getvalue(),
        )
        self.assertNotIn("style", output.getvalue())

    def test_public_build_function_uses_packaged_rotamers(self):
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            atoms = pepad_tools.build_sheets(
                "GNNQQNY",
                "package_class1.pdb",
                classes=1,
                chains=2,
            )

        self.assertIsInstance(atoms, pd.DataFrame)
        self.assertFalse(atoms.empty)
        coordinates = atoms[["x", "y", "z"]].to_numpy(dtype=float)
        self.assertTrue(np.isfinite(coordinates).all())
        self.assertTrue(Path("package_class1.pdb").is_file())
        self.assertIn("Sheets structure: Class 1", output.getvalue())


if __name__ == "__main__":
    unittest.main()
