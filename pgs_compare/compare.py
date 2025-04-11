"""
Main module for the PGS Compare package.
"""

import os
import logging

from pgs_compare.download import setup_environment
from pgs_compare.calculate import run_pgs_calculation
from pgs_compare.analyze import analyze_scores
from pgs_compare.visualize import visualize_analysis

# Set up logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


class PGSCompare:
    """
    Main class for comparing PGS scores across ancestry groups.
    """

    def __init__(self, data_dir=None, download_data=True):
        """
        Initialize the PGSCompare class.

        Args:
            data_dir (str, optional): Directory to store data. Default is "data" in the current directory.
            download_data (bool): Whether to download missing data during initialization. Default is True.
        """
        # Ensure data_dir is an absolute path if provided
        if data_dir:
            self.data_dir = os.path.abspath(data_dir)
        else:
            self.data_dir = os.path.join(os.getcwd(), "data")

        self.genomes_dir = os.path.join(self.data_dir, "1000_genomes")
        self.reference_dir = os.path.join(self.data_dir, "reference")
        self.results_dir = os.path.join(os.getcwd(), "results")

        # Ensure results_dir is an absolute path
        self.results_dir = os.path.abspath(self.results_dir)

        # Create directories
        for directory in [
            self.data_dir,
            self.genomes_dir,
            self.reference_dir,
            self.results_dir,
        ]:
            os.makedirs(directory, exist_ok=True)

        # Set up environment automatically
        logger.info("Checking dependencies and data availability")
        self.setup_results = setup_environment(
            data_dir=self.data_dir, download_data=download_data
        )

        # Log setup results
        if self.setup_results["plink_installed"]:
            logger.info("PLINK2 is installed")
        else:
            logger.warning("PLINK2 is not installed or not in PATH")

        if self.setup_results["nextflow_installed"]:
            logger.info("Nextflow is installed")
        else:
            logger.warning("Nextflow is not installed or not in PATH")

        if self.setup_results["pgsc_calc_installed"]:
            logger.info("pgsc_calc is installed/updated")

        if self.setup_results["1000_genomes_downloaded"]:
            logger.info("1000 Genomes data is available")
        else:
            logger.warning("1000 Genomes data is missing")

        if self.setup_results["reference_panels_downloaded"]:
            logger.info("Reference panels are available")
        else:
            logger.warning("Reference panels are missing")

    def calculate(
        self,
        trait_id,
        include_child_pgs=True,
        max_variants=None,
        run_ancestry=False,
        reference_panel=None,
    ):
        """
        Run PGS calculations for a specific trait.

        Args:
            trait_id (str): Trait ID (e.g., "MONDO_0005180" for Parkinson's disease)
            include_child_pgs (bool): Whether to include child-associated PGS IDs
            max_variants (int, optional): Maximum number of variants to include in PGS
            run_ancestry (bool): Whether to run ancestry analysis
            reference_panel (str, optional): Path to reference panel for ancestry analysis.
                If None and run_ancestry is True, uses the default reference panel.

        Returns:
            dict: Information about the calculation
        """
        output_dir = os.path.join(self.results_dir, trait_id)

        # Check if 1000 Genomes data is available
        panel_file = os.path.join(
            self.genomes_dir, "integrated_call_samples_v3.20130502.ALL.panel"
        )
        if not os.path.exists(panel_file):
            logger.warning(
                "1000 Genomes data not found. Run setup() or download manually."
            )

        # Run PGS calculation
        return run_pgs_calculation(
            trait_id=trait_id,
            data_dir=self.genomes_dir,
            output_dir=output_dir,
            include_child_pgs=include_child_pgs,
            max_variants=max_variants,
            run_ancestry=run_ancestry,
            reference_panel=reference_panel
            or os.path.join(self.reference_dir, "pgsc_1000G_v1.tar.zst"),
        )

    def analyze(self, trait_id=None, scores_file=None):
        """
        Analyze PGS scores across ancestry groups.

        Args:
            trait_id (str, optional): Trait ID. Used for organizing output if provided.
            scores_file (str, optional): Path to the scores file (aggregated_scores.txt.gz).
                If None, will look in the standard location based on trait_id.

        Returns:
            dict: Analysis results
        """
        # Determine output directory
        if trait_id:
            output_dir = os.path.join(self.results_dir, trait_id, "analysis")
        else:
            output_dir = os.path.join(self.results_dir, "analysis")

        return analyze_scores(
            trait_id=trait_id,
            scores_file=scores_file,
            data_dir=self.genomes_dir,
            output_dir=output_dir,
        )

    def visualize(self, trait_id=None, analysis_results=None, show_error_bars=False):
        """
        Visualize PGS analysis results.

        Args:
            trait_id (str, optional): Trait ID. Used for organizing output if provided.
            analysis_results (dict, optional): Analysis results from analyze().
                If None, will try to load from the standard location based on trait_id.
            show_error_bars (bool): Whether to show error bars for all plots. Default is False.

        Returns:
            dict: Dictionary with paths to the generated plots
        """
        # Determine directories
        if trait_id:
            analysis_dir = os.path.join(self.results_dir, trait_id, "analysis")
            output_dir = os.path.join(self.results_dir, trait_id, "analysis", "plots")
        else:
            analysis_dir = os.path.join(self.results_dir, "analysis")
            output_dir = os.path.join(self.results_dir, "analysis", "plots")

        return visualize_analysis(
            analysis_results=analysis_results,
            analysis_dir=analysis_dir,
            output_dir=output_dir,
            show_error_bars=show_error_bars,
        )

    def run_pipeline(
        self,
        trait_id,
        include_child_pgs=True,
        max_variants=None,
        run_ancestry=False,
        visualize=True,
        show_error_bars=False,
    ):
        """
        Run the full pipeline (calculate, analyze, visualize) for a specific trait.

        Args:
            trait_id (str): Trait ID (e.g., "MONDO_0005180" for Parkinson's disease)
            include_child_pgs (bool): Whether to include child-associated PGS IDs
            max_variants (int, optional): Maximum number of variants to include in PGS
            run_ancestry (bool): Whether to run ancestry analysis
            visualize (bool): Whether to generate visualization plots
            show_error_bars (bool): Whether to show error bars for all plots. Default is False.

        Returns:
            dict: Pipeline results
        """
        # Run calculation
        calc_results = self.calculate(
            trait_id=trait_id,
            include_child_pgs=include_child_pgs,
            max_variants=max_variants,
            run_ancestry=run_ancestry,
        )

        if not calc_results["success"]:
            return {
                "success": False,
                "stage": "calculation",
                "error": "PGS calculation failed",
                "calculation_results": calc_results,
            }

        # Run analysis
        analysis_results = self.analyze(trait_id=trait_id)

        if not analysis_results["success"]:
            return {
                "success": False,
                "stage": "analysis",
                "error": "Analysis failed",
                "calculation_results": calc_results,
                "analysis_results": analysis_results,
            }

        # Run visualization if requested
        if visualize:
            viz_results = self.visualize(
                trait_id=trait_id,
                analysis_results=analysis_results,
                show_error_bars=show_error_bars,
            )

            return {
                "success": True,
                "calculation_results": calc_results,
                "analysis_results": analysis_results,
                "visualization_results": viz_results,
            }

        return {
            "success": True,
            "calculation_results": calc_results,
            "analysis_results": analysis_results,
        }
