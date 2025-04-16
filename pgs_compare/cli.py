"""
Command-line interface for PGS Compare.
"""

import argparse
import logging
import sys
import os
from pgs_compare.compare import PGSCompare

logger = logging.getLogger(__name__)


def setup_parser():
    """
    Set up the command-line argument parser.

    Returns:
        argparse.ArgumentParser: Configured argument parser
    """
    parser = argparse.ArgumentParser(
        description="PGS Compare - Compare Polygenic Scores across ancestry groups"
    )

    # Common optional arguments
    parser.add_argument(
        "--data-dir", help="Directory for data storage (default: ./data)"
    )
    parser.add_argument(
        "--results-dir", help="Directory for results (default: ./results)"
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Enable verbose output"
    )

    # Sub-parsers for commands
    subparsers = parser.add_subparsers(dest="command", help="Command to run")

    # Calculate command
    calc_parser = subparsers.add_parser(
        "calculate", help="Run PGS calculations for a specific trait"
    )
    calc_parser.add_argument(
        "trait_id", help="Trait ID (e.g., MONDO_0005180 for Parkinson's disease)"
    )
    calc_parser.add_argument(
        "--exclude-child-pgs",
        action="store_true",
        help="Exclude child-associated PGS IDs",
    )
    calc_parser.add_argument(
        "--max-variants",
        type=int,
        default=1500000,
        help="Maximum number of variants to include in PGS (default: 1500000)",
    )
    calc_parser.add_argument(
        "--run-ancestry", action="store_true", help="Run ancestry analysis"
    )
    calc_parser.add_argument(
        "--reference-panel", help="Path to reference panel for ancestry analysis"
    )
    calc_parser.add_argument(
        "--pgs-ids",
        help="Custom comma-separated string of PGS IDs to calculate (e.g., PGS001229,PGS001405). If provided, trait_id is optional.",
    )

    # Analyze command
    analyze_parser = subparsers.add_parser(
        "analyze", help="Analyze PGS scores across ancestry groups"
    )
    analyze_parser.add_argument("--trait-id", help="Trait ID for organizing output")
    analyze_parser.add_argument(
        "--scores-file", help="Path to the scores file (aggregated_scores.txt.gz)"
    )

    # Visualize command
    visualize_parser = subparsers.add_parser(
        "visualize", help="Visualize PGS analysis results"
    )
    visualize_parser.add_argument("--trait-id", help="Trait ID for organizing output")
    visualize_parser.add_argument(
        "--show-error-bars",
        action="store_true",
        help="Show error bars on plots (default: False)",
    )

    # Pipeline command
    pipeline_parser = subparsers.add_parser(
        "pipeline", help="Run the full pipeline (calculate, analyze, visualize)"
    )
    pipeline_parser.add_argument(
        "trait_id", help="Trait ID (e.g., MONDO_0005180 for Parkinson's disease)"
    )
    pipeline_parser.add_argument(
        "--data-dir", help="Directory for data storage (default: ./data)"
    )
    pipeline_parser.add_argument(
        "--exclude-child-pgs",
        action="store_true",
        help="Exclude child-associated PGS IDs",
    )
    pipeline_parser.add_argument(
        "--max-variants",
        type=int,
        default=1500000,
        help="Maximum number of variants to include in PGS (default: 1500000)",
    )
    pipeline_parser.add_argument(
        "--run-ancestry", action="store_true", help="Run ancestry analysis"
    )
    pipeline_parser.add_argument(
        "--skip-visualize", action="store_true", help="Skip visualization step"
    )
    pipeline_parser.add_argument(
        "--show-error-bars",
        action="store_true",
        help="Show error bars on plots (default: False)",
    )
    pipeline_parser.add_argument(
        "--pgs-ids",
        help="Custom comma-separated string of PGS IDs to calculate (e.g., PGS001229,PGS001405). If provided, trait_id is optional.",
    )

    return parser


def main():
    """
    Main entry point for the command-line interface.

    Returns:
        int: Exit code (0 for success, non-zero for errors)
    """
    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    # Parse command-line arguments
    parser = setup_parser()
    args = parser.parse_args()

    # Display help if no command is specified
    if args.command is None:
        parser.print_help()
        return 0

    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Initialize PGSCompare
    pgs_compare = PGSCompare(
        data_dir=args.data_dir or os.path.join(os.getcwd(), "data")
    )

    # Override results_dir if provided
    if hasattr(args, "results_dir") and args.results_dir:
        pgs_compare.results_dir = args.results_dir

    # Check if required dependencies are available
    if (
        not pgs_compare.setup_results["plink_installed"]
        or not pgs_compare.setup_results["nextflow_installed"]
    ):
        logger.error(
            "Some required dependencies are missing. Please install them and try again."
        )
        return 1

    # Execute the requested command
    if args.command == "calculate":
        result = pgs_compare.calculate(
            trait_id=args.trait_id,
            include_child_pgs=not args.exclude_child_pgs,
            max_variants=args.max_variants,
            run_ancestry=args.run_ancestry,
            reference_panel=args.reference_panel,
            pgs_ids=args.pgs_ids,
        )

        if not result["success"]:
            logger.error("Calculation failed")
            return 1

        logger.info(f"Calculation completed for {len(result['pgs_ids'])} PGS IDs")
        logger.info(f"Results saved to: {result['output_path']}")

    elif args.command == "analyze":
        result = pgs_compare.analyze(
            trait_id=args.trait_id, scores_file=args.scores_file
        )

        if not result["success"]:
            logger.error(f"Analysis failed: {result.get('error', 'Unknown error')}")
            return 1

        logger.info(f"Analysis completed")
        logger.info(f"Results saved to: {result['output_path']}")

    elif args.command == "visualize":
        result = pgs_compare.visualize(
            trait_id=args.trait_id, show_error_bars=args.show_error_bars
        )

        if not result["success"]:
            logger.error(
                f"Visualization failed: {result.get('error', 'Unknown error')}"
            )
            return 1

        logger.info(f"Visualization completed")
        logger.info(
            f"Plots saved to: {os.path.dirname(next(iter(result['plots'].values()), ''))}"
        )

    elif args.command == "pipeline":
        result = pgs_compare.run_pipeline(
            trait_id=args.trait_id,
            include_child_pgs=not args.exclude_child_pgs,
            max_variants=args.max_variants,
            run_ancestry=args.run_ancestry,
            visualize=not args.skip_visualize,
            show_error_bars=args.show_error_bars,
            pgs_ids=args.pgs_ids,
        )

        if not result["success"]:
            logger.error(f"Pipeline failed at stage: {result.get('stage', 'unknown')}")
            return 1

        logger.info(f"Pipeline completed successfully")

        # Print paths to results
        if "calculation_results" in result:
            logger.info(
                f"Calculation results: {result['calculation_results']['output_path']}"
            )
        if "analysis_results" in result:
            logger.info(
                f"Analysis results: {result['analysis_results']['output_path']}"
            )

    return 0


if __name__ == "__main__":
    sys.exit(main())
