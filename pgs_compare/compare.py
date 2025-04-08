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

    def __init__(self, data_dir=None, download_data=False):
        """
        Initialize the PGSCompare class.

        Args:
            data_dir (str, optional): Directory to store data. Default is "data" in the current directory.
            download_data (bool): Whether to download data during initialization.
        """
        self.data_dir = data_dir or os.path.join(os.getcwd(), "data")
        self.genomes_dir = os.path.join(self.data_dir, "1000_genomes")
        self.reference_dir = os.path.join(self.data_dir, "reference")
        self.results_dir = os.path.join(os.getcwd(), "results")

        # Create directories
        for directory in [
            self.data_dir,
            self.genomes_dir,
            self.reference_dir,
            self.results_dir,
        ]:
            os.makedirs(directory, exist_ok=True)

        # Set up environment if requested
        if download_data:
            self.setup()

    def setup(self, download_genomes=True, download_reference=True):
        """
        Set up the environment for PGS comparisons.

        Args:
            download_genomes (bool): Whether to download 1000 Genomes data.
            download_reference (bool): Whether to download reference panels.

        Returns:
            dict: Status of each setup step
        """
        logger.info("Setting up PGS Compare environment")

        setup_results = setup_environment(
            data_dir=self.data_dir,
            download_data=False,  # We'll handle downloading separately
        )

        # Download 1000 Genomes data if requested
        if download_genomes and setup_results["plink_installed"]:
            from pgs_compare.download import download_1000_genomes

            setup_results["1000_genomes_downloaded"] = download_1000_genomes(
                data_dir=self.genomes_dir
            )

        # Download reference panels if requested
        if download_reference:
            from pgs_compare.download import download_reference_panels

            setup_results["reference_panels_downloaded"] = download_reference_panels(
                data_dir=self.reference_dir
            )

        return setup_results

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

    def visualize(self, trait_id=None, analysis_results=None):
        """
        Visualize PGS analysis results.

        Args:
            trait_id (str, optional): Trait ID. Used for organizing output if provided.
            analysis_results (dict, optional): Analysis results from analyze().
                If None, will try to load from the standard location based on trait_id.

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
        )

    def run_pipeline(
        self,
        trait_id,
        include_child_pgs=True,
        max_variants=None,
        run_ancestry=False,
        visualize=True,
    ):
        """
        Run the full pipeline (calculate, analyze, visualize) for a specific trait.

        Args:
            trait_id (str): Trait ID (e.g., "MONDO_0005180" for Parkinson's disease)
            include_child_pgs (bool): Whether to include child-associated PGS IDs
            max_variants (int, optional): Maximum number of variants to include in PGS
            run_ancestry (bool): Whether to run ancestry analysis
            visualize (bool): Whether to generate visualization plots

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
                trait_id=trait_id, analysis_results=analysis_results
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
