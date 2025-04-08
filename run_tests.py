#!/usr/bin/env python
"""
Script to run tests for the pgs_compare package.
"""

import unittest
import sys

if __name__ == "__main__":
    test_loader = unittest.TestLoader()
    test_suite = test_loader.discover("tests", pattern="test_*.py")
    test_runner = unittest.TextTestRunner(verbosity=2)
    result = test_runner.run(test_suite)
    sys.exit(not result.wasSuccessful())
