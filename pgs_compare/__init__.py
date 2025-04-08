"""
PGS Compare - A package for comparing Polygenic Scores across ancestry groups.
"""

from pgs_compare.compare import PGSCompare
from pgs_compare.calculate import run_pgs_calculation
from pgs_compare.analyze import analyze_scores

__version__ = "0.1.0"
__all__ = ["PGSCompare", "run_pgs_calculation", "analyze_scores"]
