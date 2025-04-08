"""
Module for downloading necessary data for PGS comparisons.
"""

import os
import subprocess
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def check_plink_installation():
    """
    Check if PLINK2 is installed and accessible.
    Returns True if PLINK2 is installed, False otherwise.
    """
    try:
        result = subprocess.run(
            ["plink2", "--version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        if result.returncode == 0:
            logger.info(f"PLINK2 is installed: {result.stdout.strip()}")
            return True
        else:
            logger.warning("PLINK2 check failed with non-zero return code")
            return False
    except FileNotFoundError:
        logger.warning("PLINK2 is not installed or not in PATH")
        return False


def check_nextflow_installation():
    """
    Check if Nextflow is installed and accessible.
    Returns True if Nextflow is installed, False otherwise.
    """
    try:
        result = subprocess.run(
            ["nextflow", "-version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        if result.returncode == 0:
            logger.info(f"Nextflow is installed: {result.stdout.strip()}")
            return True
        else:
            logger.warning("Nextflow check failed with non-zero return code")
            return False
    except FileNotFoundError:
        logger.warning("Nextflow is not installed or not in PATH")
        return False


def download_1000_genomes(data_dir=None):
    """
    Download 1000 Genomes data.

    Args:
        data_dir (str, optional): Path to store the data.
                                 Default is 'data/1000_genomes' in the current directory.

    Returns:
        bool: True if download is successful, False otherwise.
    """
    if not check_plink_installation():
        logger.error(
            "PLINK2 is required for 1000 Genomes data processing but is not installed."
        )
        return False

    if data_dir is None:
        data_dir = os.path.join(os.getcwd(), "data", "1000_genomes")

    script_path = os.path.join(
        os.path.dirname(__file__), "scripts", "download_1000_genomes.sh"
    )

    # Make script executable
    os.chmod(script_path, 0o755)

    try:
        logger.info(f"Starting 1000 Genomes download to {data_dir}")
        subprocess.run([script_path, data_dir], check=True)
        logger.info("1000 Genomes download completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Error downloading 1000 Genomes data: {e}")
        return False


def download_reference_panels(data_dir=None):
    """
    Download reference panels for ancestry analysis.

    Args:
        data_dir (str, optional): Path to store the reference panels.
                                 Default is 'data/reference' in the current directory.

    Returns:
        bool: True if download is successful, False otherwise.
    """
    if data_dir is None:
        data_dir = os.path.join(os.getcwd(), "data", "reference")

    script_path = os.path.join(
        os.path.dirname(__file__), "scripts", "download_reference.sh"
    )

    # Make script executable
    os.chmod(script_path, 0o755)

    try:
        logger.info(f"Starting reference panels download to {data_dir}")
        subprocess.run([script_path, data_dir], check=True)
        logger.info("Reference panels download completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Error downloading reference panels: {e}")
        return False


def install_pgsc_calc():
    """
    Install or pull the latest version of pgsc_calc.

    Returns:
        bool: True if installation/update is successful, False otherwise.
    """
    if not check_nextflow_installation():
        logger.error("Nextflow is required for pgsc_calc but is not installed.")
        return False

    try:
        logger.info("Pulling latest pgsc_calc")
        subprocess.run(["nextflow", "pull", "pgscatalog/pgsc_calc"], check=True)
        logger.info("pgsc_calc successfully updated")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Error updating pgsc_calc: {e}")
        return False


def setup_environment(data_dir=None, download_data=True):
    """
    Set up the environment for PGS comparisons.

    Args:
        data_dir (str, optional): Base directory for data storage.
        download_data (bool): Whether to download data or not.

    Returns:
        dict: Status of each setup step
    """
    results = {
        "plink_installed": check_plink_installation(),
        "nextflow_installed": check_nextflow_installation(),
        "pgsc_calc_installed": False,
        "1000_genomes_downloaded": False,
        "reference_panels_downloaded": False,
    }

    # Create data directory structure
    if data_dir is None:
        data_dir = os.path.join(os.getcwd(), "data")

    os.makedirs(data_dir, exist_ok=True)
    genomes_dir = os.path.join(data_dir, "1000_genomes")
    reference_dir = os.path.join(data_dir, "reference")
    
    # Create subdirectories
    os.makedirs(genomes_dir, exist_ok=True)
    os.makedirs(reference_dir, exist_ok=True)

    # Install pgsc_calc
    if results["nextflow_installed"]:
        results["pgsc_calc_installed"] = install_pgsc_calc()

    # Check if 1000 Genomes data exists
    panel_file = os.path.join(genomes_dir, "integrated_call_samples_v3.20130502.ALL.panel")
    has_genomes_data = os.path.exists(panel_file)
    
    # Check if reference panels exist
    reference_file = os.path.join(reference_dir, "pgsc_1000G_v1.tar.zst")
    has_reference_data = os.path.exists(reference_file)
    
    # Download data if requested and not already present
    if download_data:
        if not has_genomes_data and results["plink_installed"]:
            results["1000_genomes_downloaded"] = download_1000_genomes(genomes_dir)
        else:
            results["1000_genomes_downloaded"] = has_genomes_data
            if has_genomes_data:
                logger.info("1000 Genomes data already exists, skipping download")

        if not has_reference_data:
            results["reference_panels_downloaded"] = download_reference_panels(reference_dir)
        else:
            results["reference_panels_downloaded"] = has_reference_data
            if has_reference_data:
                logger.info("Reference panels already exist, skipping download")

    return results
