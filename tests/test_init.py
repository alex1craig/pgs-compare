"""
Tests for package initialization.
"""

import os
import tempfile
import unittest
from pgs_compare import PGSCompare


class TestPGSCompare(unittest.TestCase):
    """
    Tests for the PGSCompare class.
    """

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.TemporaryDirectory()
        self.data_dir = os.path.join(self.temp_dir.name, "data")

    def tearDown(self):
        """Tear down test fixtures."""
        self.temp_dir.cleanup()

    def test_init(self):
        """Test PGSCompare initialization."""
        pgs_compare = PGSCompare(data_dir=self.data_dir)

        # Check that directories were created
        self.assertTrue(os.path.exists(pgs_compare.data_dir))
        self.assertTrue(os.path.exists(pgs_compare.genomes_dir))
        self.assertTrue(os.path.exists(pgs_compare.reference_dir))
        self.assertTrue(os.path.exists(pgs_compare.results_dir))

        # Check directory structure
        self.assertEqual(pgs_compare.data_dir, self.data_dir)
        self.assertEqual(
            pgs_compare.genomes_dir, os.path.join(self.data_dir, "1000_genomes")
        )
        self.assertEqual(
            pgs_compare.reference_dir, os.path.join(self.data_dir, "reference")
        )


if __name__ == "__main__":
    unittest.main()
