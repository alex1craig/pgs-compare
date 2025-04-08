"""
Module for running PGS calculations on 1000 Genomes data.
"""

import os
import subprocess
import logging
import requests
import json
import pandas as pd
import shutil
from pathlib import Path

logger = logging.getLogger(__name__)


def fetch_pgs_ids(trait_id, include_child_pgs=True, max_variants=None):
    """
    Fetch PGS IDs for a specific trait from the PGS Catalog API.

    Args:
        trait_id (str): Trait ID (e.g., "MONDO_0005180" for Parkinson's disease)
        include_child_pgs (bool): Whether to include child-associated PGS IDs
        max_variants (int, optional): Maximum number of variants for PGS to include.
                                     If None, all PGS scores are included.

    Returns:
        tuple: (list of PGS IDs, trait label, dict of trait response)
    """
    logger.info(f"Fetching PGS IDs for trait {trait_id}")

    # Make API request to get trait information
    try:
        trait_response = requests.get(
            f"https://www.pgscatalog.org/rest/trait/{trait_id}?include_children=1"
        ).json()
    except requests.RequestException as e:
        logger.error(f"Error fetching trait information: {e}")
        return [], None, {}

    # Set up PGS IDs list from the response
    pgs_ids = []
    if "associated_pgs_ids" in trait_response:
        pgs_ids.extend(trait_response["associated_pgs_ids"])

    # Include child-associated PGS IDs if requested
    if include_child_pgs and "child_associated_pgs_ids" in trait_response:
        pgs_ids.extend(trait_response["child_associated_pgs_ids"])
        # For compatibility with API structure
        trait_response["associated_pgs_ids"] = pgs_ids

    trait_label = trait_response.get("label", "")

    # Filter by variant count if max_variants is specified
    if max_variants is not None:
        filtered_pgs_ids = []
        for pgs_id in pgs_ids:
            try:
                score_response = requests.get(
                    f"https://www.pgscatalog.org/rest/score/{pgs_id}"
                ).json()
                variant_count = score_response.get("variants_number", 0)

                if variant_count <= max_variants:
                    filtered_pgs_ids.append(pgs_id)
                else:
                    logger.info(
                        f"Skipping {pgs_id} due to high variant count: {variant_count} > {max_variants}"
                    )
            except requests.RequestException as e:
                logger.warning(f"Error fetching data for {pgs_id}: {e}")
                continue

        pgs_ids = filtered_pgs_ids

    logger.info(f"Found {len(pgs_ids)} PGS IDs for trait {trait_id} ({trait_label})")
    return pgs_ids, trait_label, trait_response


def create_samplesheet(data_dir, chromosomes=None, output_dir=None):
    """
    Create a samplesheet for PGS calculation.

    Args:
        data_dir (str): Directory containing the 1000 Genomes data
        chromosomes (list, optional): List of chromosomes to include. Default is 1-22.
        output_dir (str, optional): Directory to save the samplesheet. Default is data_dir/pgs.

    Returns:
        str: Path to the generated samplesheet
    """
    if chromosomes is None:
        chromosomes = [str(i) for i in range(1, 23)]

    if output_dir is None:
        output_dir = os.path.join(data_dir, "pgs")

    os.makedirs(output_dir, exist_ok=True)

    # Create DataFrame for samplesheet
    samplesheet_data = []
    for chromosome in chromosomes:
        samplesheet_data.append(
            {
                "sampleset": "ALL",
                "path_prefix": os.path.join(data_dir, f"chr_{chromosome}"),
                "chrom": chromosome,
                "format": "pfile",
            }
        )

    samplesheet = pd.DataFrame(samplesheet_data)
    samplesheet_path = os.path.join(output_dir, "samplesheet.csv")
    samplesheet.to_csv(samplesheet_path, index=False)

    logger.info(f"Created samplesheet at {samplesheet_path}")
    return samplesheet_path


def run_pgs_calculation(
    trait_id,
    data_dir=None,
    output_dir=None,
    include_child_pgs=True,
    max_variants=1500000,
    run_ancestry=False,
    reference_panel=None,
):
    """
    Run PGS calculations for a specific trait.

    Args:
        trait_id (str): Trait ID (e.g., "MONDO_0005180" for Parkinson's disease)
        data_dir (str, optional): Directory containing the 1000 Genomes data
        output_dir (str, optional): Directory to store the output
        include_child_pgs (bool): Whether to include child-associated PGS IDs
        max_variants (int, optional): Maximum number of variants to include in PGS
        run_ancestry (bool): Whether to run ancestry analysis
        reference_panel (str, optional): Path to reference panel for ancestry analysis.
            If None and run_ancestry is True, uses the default reference panel.

    Returns:
        dict: Information about the calculation, including:
            - pgs_ids: List of PGS IDs included in the calculation
            - trait_label: Label for the trait
            - output_path: Path to the output directory
            - success: Whether the calculation was successful
    """
    # Set up directories
    if data_dir is None:
        data_dir = os.path.join(os.getcwd(), "data", "1000_genomes")

    if output_dir is None:
        output_dir = os.path.join(os.getcwd(), "results")

    os.makedirs(output_dir, exist_ok=True)

    # Fetch PGS IDs for the trait
    pgs_ids, trait_label, trait_response = fetch_pgs_ids(
        trait_id, include_child_pgs=include_child_pgs, max_variants=max_variants
    )

    if not pgs_ids:
        logger.error(f"No valid PGS IDs found for trait {trait_id}")
        return {
            "pgs_ids": [],
            "trait_label": trait_label,
            "output_path": output_dir,
            "success": False,
        }

    # Create samplesheet
    samplesheet_path = create_samplesheet(
        data_dir, output_dir=os.path.join(output_dir, trait_id, "input")
    )

    # Prepare Nextflow command
    cmd = [
        "nextflow",
        "run",
        "pgscatalog/pgsc_calc",
        "-profile",
        "docker",
        "--input",
        samplesheet_path,
        "--target_build",
        "GRCh37",
        "--pgs_id",
        ",".join(pgs_ids),
    ]

    # Add ancestry analysis if requested
    if run_ancestry:
        if reference_panel is None:
            reference_panel = os.path.join(
                os.path.dirname(data_dir), "reference", "pgsc_1000G_v1.tar.zst"
            )

        if os.path.exists(reference_panel):
            cmd.extend(["--run_ancestry", reference_panel])
            logger.info(f"Using reference panel: {reference_panel}")
        else:
            logger.warning(f"Reference panel not found: {reference_panel}")
            logger.warning("Proceeding without ancestry analysis")

    # Run the command
    logger.info(f"Running PGS calculation for trait {trait_id} ({trait_label})")
    logger.info(f"Command: {' '.join(cmd)}")

    try:
        # Set the working directory to ensure output goes to the right place
        process = subprocess.run(
            cmd,
            check=True,
            cwd=(
                os.path.join(output_dir, trait_id)
                if os.path.exists(os.path.join(output_dir, trait_id))
                else output_dir
            ),
        )

        # Copy the results to a more accessible location
        result_dir = os.path.join(output_dir, trait_id, "results")
        os.makedirs(result_dir, exist_ok=True)

        # Determine source path based on where nextflow put the results
        # If in output_dir/ALL/score
        all_score_dir = os.path.join(output_dir, trait_id, "ALL", "score")
        results_score_dir = os.path.join(output_dir, trait_id, "results", "score")

        if os.path.exists(all_score_dir):
            # Copy from ALL/score to results/score
            if os.path.exists(results_score_dir):
                shutil.rmtree(results_score_dir)
            shutil.copytree(all_score_dir, results_score_dir)
        else:
            # Check for results in output_dir/results/score
            source_score_dir = os.path.join(output_dir, "results", "score")
            if os.path.exists(source_score_dir):
                if os.path.exists(results_score_dir):
                    shutil.rmtree(results_score_dir)
                shutil.copytree(source_score_dir, results_score_dir)
            else:
                logger.warning(f"Could not find score results to copy")

        # Save trait information
        trait_info_path = os.path.join(result_dir, f"{trait_id}_info.json")
        with open(trait_info_path, "w") as f:
            json.dump(trait_response, f)

        logger.info(f"PGS calculation completed for trait {trait_id}")
        return {
            "pgs_ids": pgs_ids,
            "trait_label": trait_label,
            "output_path": (
                results_score_dir if os.path.exists(results_score_dir) else result_dir
            ),
            "success": True,
        }

    except subprocess.CalledProcessError as e:
        logger.error(f"Error running PGS calculation: {e}")
        return {
            "pgs_ids": pgs_ids,
            "trait_label": trait_label,
            "output_path": output_dir,
            "success": False,
        }
