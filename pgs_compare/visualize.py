"""
Module for visualizing PGS analysis results.
"""

import os
import logging
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

logger = logging.getLogger(__name__)


def setup_plot(title, xlabel, ylabel, add_zero_line=False, trait_name=None):
    """
    Set up plot with common formatting.

    Args:
        title (str): Plot title
        xlabel (str): X-axis label
        ylabel (str): Y-axis label
        add_zero_line (bool): Whether to add a horizontal line at y=0
        trait_name (str, optional): Trait name to add to title

    Returns:
        None
    """
    if trait_name:
        title = f"{title} ({trait_name})"

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()

    if add_zero_line:
        plt.axvline(x=0, linestyle="--", color="black", alpha=0.5)
        plt.grid(alpha=0.3)

    plt.tight_layout()


def plot_kde(scores, label, alpha=0.75):
    """
    Plot kernel density estimate for the given scores.

    Args:
        scores (pandas.Series): Scores to plot
        label (str): Label for the plot
        alpha (float): Alpha value for the plot

    Returns:
        bool: True if the plot was created, False otherwise
    """
    if len(scores) > 0:
        density = stats.gaussian_kde(scores)
        xs = np.linspace(scores.min(), scores.max(), 200)
        plt.plot(xs, density(xs), alpha=alpha, label=label)
        return True
    return False


def plot_distribution_by_ancestry(scores_df, pgs_id, output_dir=None, trait_name=None):
    """
    Plot distribution of scores by ancestry group for a specific PGS.

    Args:
        scores_df (pandas.DataFrame): DataFrame with scores and ancestry information
        pgs_id (str): PGS ID to plot
        output_dir (str, optional): Directory to save the plot
        trait_name (str, optional): Trait name to add to title

    Returns:
        str or None: Path to the saved plot, or None if no plot was created
    """
    plt.figure(figsize=(10, 6))

    # Get unique ancestry groups
    ancestry_groups = list(scores_df["GROUP"].unique()) + ["ALL"]

    # Plot distribution for each ancestry group
    for group in ancestry_groups:
        # Get scores for this ancestry group
        group_scores = scores_df[scores_df["PGS"] == pgs_id]
        scores = (
            group_scores["SUM"]
            if group == "ALL"
            else group_scores[group_scores["GROUP"] == group]["SUM"]
        )

        # Plot distribution
        plot_kde(scores, group, alpha=1)

    # Set up plot
    setup_plot(
        f"Distribution of PGS Scores by Ancestry for {pgs_id}",
        "Polygenic Score",
        "Density",
        trait_name=trait_name,
    )

    # Save plot
    if output_dir:
        os.makedirs(os.path.join(output_dir, "distributions"), exist_ok=True)
        plot_path = os.path.join(
            output_dir, "distributions", f"{pgs_id}_distributions.png"
        )
        plt.savefig(plot_path)
        plt.close()
        return plot_path
    else:
        plt.show()
        plt.close()
        return None


def plot_distribution_by_pgs(scores_df, group, output_dir=None, trait_name=None):
    """
    Plot distribution of scores by PGS for a specific ancestry group.

    Args:
        scores_df (pandas.DataFrame): DataFrame with scores and ancestry information
        group (str): Ancestry group to plot
        output_dir (str, optional): Directory to save the plot
        trait_name (str, optional): Trait name to add to title

    Returns:
        str or None: Path to the saved plot, or None if no plot was created
    """
    plt.figure(figsize=(10, 6))

    # Get scores for this ancestry group
    group_scores = (
        scores_df if group == "ALL" else scores_df[scores_df["GROUP"] == group]
    )

    # Plot distribution for each PGS
    for pgs_id in scores_df["PGS"].unique():
        # Get scores for this PGS
        scores = group_scores[group_scores["PGS"] == pgs_id]["SUM"]

        # Plot distribution
        plot_kde(scores, pgs_id, alpha=1)

    # Set up plot
    setup_plot(
        f"Distribution of PGS Scores in {group}",
        "Polygenic Score",
        "Density",
        trait_name=trait_name,
    )

    # Save plot
    if output_dir:
        os.makedirs(os.path.join(output_dir, "distributions"), exist_ok=True)
        plot_path = os.path.join(
            output_dir, "distributions", f"{group}_distributions.png"
        )
        plt.savefig(plot_path)
        plt.close()
        return plot_path
    else:
        plt.show()
        plt.close()
        return None


def plot_standardized_distribution_by_ancestry(
    scores_df, pgs_id, output_dir=None, trait_name=None
):
    """
    Plot standardized distribution of scores by ancestry group for a specific PGS.

    Args:
        scores_df (pandas.DataFrame): DataFrame with scores and ancestry information
        pgs_id (str): PGS ID to plot
        output_dir (str, optional): Directory to save the plot
        trait_name (str, optional): Trait name to add to title

    Returns:
        str or None: Path to the saved plot, or None if no plot was created
    """
    plt.figure(figsize=(12, 6))

    # Get unique ancestry groups
    ancestry_groups = list(scores_df["GROUP"].unique()) + ["ALL"]

    # Plot distribution for each ancestry group
    for group in ancestry_groups:
        # Get z-scores for this group
        mask = (
            scores_df["PGS"] == pgs_id
            if group == "ALL"
            else ((scores_df["PGS"] == pgs_id) & (scores_df["GROUP"] == group))
        )
        group_z_scores = scores_df.loc[mask, "z_score"]

        # Plot distribution
        plot_kde(group_z_scores, group, alpha=1)

    # Set up plot
    setup_plot(
        f"Standardized PGS Scores for {pgs_id}",
        "Z-Score",
        "Density",
        add_zero_line=True,
        trait_name=trait_name,
    )

    # Save plot
    if output_dir:
        os.makedirs(
            os.path.join(output_dir, "standardized_distributions"), exist_ok=True
        )
        plot_path = os.path.join(
            output_dir,
            "standardized_distributions",
            f"{pgs_id}_standardized_distributions.png",
        )
        plt.savefig(plot_path)
        plt.close()
        return plot_path
    else:
        plt.show()
        plt.close()
        return None


def plot_correlation_matrix(scores_df, group, output_dir=None, trait_name=None):
    """
    Plot correlation matrix for a specific ancestry group.

    Args:
        scores_df (pandas.DataFrame): DataFrame with scores and ancestry information
        group (str): Ancestry group to plot
        output_dir (str, optional): Directory to save the plot
        trait_name (str, optional): Trait name to add to title

    Returns:
        str or None: Path to the saved plot, or None if no plot was created
    """
    plt.figure(figsize=(12, 10))

    # Get scores for this ancestry group
    group_data = scores_df if group == "ALL" else scores_df[scores_df["GROUP"] == group]

    if len(group_data) == 0:
        logger.warning(f"No data for ancestry group {group}")
        plt.close()
        return None

    # Create a pivot table: rows=individuals, columns=PGS IDs, values=scores
    try:
        pivot_data = group_data.pivot(index="IID", columns="PGS", values="SUM")

        # Calculate correlation matrix
        corr_matrix = pivot_data.corr()

        # Plot heatmap
        cmap = plt.cm.RdBu_r
        im = plt.imshow(corr_matrix, cmap=cmap, vmin=-1, vmax=1)

        # Add colorbar
        cbar = plt.colorbar(im)
        cbar.set_label("Pearson Correlation")

        # Add labels and ticks
        title = f"Correlation Matrix of PGS Scores for {group} Ancestry"
        if trait_name:
            title += f" ({trait_name})"

        plt.title(title)
        plt.xticks(
            range(len(corr_matrix.columns)),
            corr_matrix.columns,
            rotation=45,
            ha="right",
        )
        plt.yticks(range(len(corr_matrix.columns)), corr_matrix.columns)

        # Add correlation values to the heatmap
        for i in range(len(corr_matrix.columns)):
            for j in range(len(corr_matrix.columns)):
                text = plt.text(
                    j,
                    i,
                    f"{corr_matrix.iloc[i, j]:.2f}",
                    ha="center",
                    va="center",
                    color="black" if abs(corr_matrix.iloc[i, j]) < 0.7 else "white",
                )

        plt.tight_layout()

        # Save plot
        if output_dir:
            os.makedirs(os.path.join(output_dir, "correlations"), exist_ok=True)
            plot_path = os.path.join(
                output_dir, "correlations", f"{group}_correlation_matrix.png"
            )
            plt.savefig(plot_path)
            plt.close()
            return plot_path
        else:
            plt.show()
            plt.close()
            return None

    except Exception as e:
        logger.warning(f"Error calculating correlation matrix for group {group}: {e}")
        plt.close()
        return None


def plot_average_correlations(average_correlations, output_dir=None, trait_name=None):
    """
    Plot average correlations across ancestry groups.

    Args:
        average_correlations (dict): Dictionary with average correlations by ancestry group
        output_dir (str, optional): Directory to save the plot
        trait_name (str, optional): Trait name to add to title

    Returns:
        str or None: Path to the saved plot, or None if no plot was created
    """
    plt.figure(figsize=(10, 6))

    # Plot average correlations
    plt.bar(
        list(average_correlations.keys()),
        [float(val) for val in average_correlations.values()],
    )

    # Set up plot
    plt.xlabel("Ancestry Group")
    plt.ylabel("Average Correlation")

    title = "Average Correlation of PGS Scores by Ancestry Group"
    if trait_name:
        title += f" ({trait_name})"

    plt.title(title)
    plt.ylim(bottom=0)
    plt.tight_layout()

    # Save plot
    if output_dir:
        os.makedirs(os.path.join(output_dir, "correlations"), exist_ok=True)
        plot_path = os.path.join(
            output_dir, "correlations", "average_correlation_bar_chart.png"
        )
        plt.savefig(plot_path)
        plt.close()
        return plot_path
    else:
        plt.show()
        plt.close()
        return None


def plot_deviations(deviations, ancestry_groups, output_dir=None, trait_name=None):
    """
    Plot mean deviations for each PGS study by ancestry group.

    Args:
        deviations (dict): Dictionary with deviation information
        ancestry_groups (list): List of ancestry groups
        output_dir (str, optional): Directory to save the plot
        trait_name (str, optional): Trait name to add to title

    Returns:
        str or None: Path to the saved plot, or None if no plot was created
    """
    plt.figure(figsize=(12, 8))

    pgs_ids = list(deviations.keys())
    x = np.arange(len(pgs_ids))
    width = 0.8 / len(ancestry_groups)

    for i, group in enumerate(ancestry_groups):
        means = [
            (
                deviations[pgs_id][group]["mean_deviation"]
                if group in deviations[pgs_id]
                else 0
            )
            for pgs_id in pgs_ids
        ]
        plt.bar(
            x + i * width - (len(ancestry_groups) - 1) * width / 2,
            means,
            width,
            label=group,
        )

    # Set up plot
    plt.xlabel("PGS Study")
    plt.ylabel("Mean Deviation from Average Z-Score")

    title = "Mean Deviation of PGS Z-Scores from Average"
    if trait_name:
        title += f" ({trait_name})"

    plt.title(title)
    plt.xticks(x, pgs_ids, rotation=45, ha="right")
    plt.legend()
    plt.axhline(y=0, color="black", linestyle="-", alpha=0.3)
    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()

    # Save plot
    if output_dir:
        os.makedirs(os.path.join(output_dir, "deviations"), exist_ok=True)
        plot_path = os.path.join(
            output_dir, "deviations", "mean_z_score_deviations_from_average.png"
        )
        plt.savefig(plot_path)
        plt.close()
        return plot_path
    else:
        plt.show()
        plt.close()
        return None


def visualize_analysis(analysis_results=None, analysis_dir=None, output_dir=None):
    """
    Visualize PGS analysis results.

    Args:
        analysis_results (dict, optional): Analysis results from analyze_scores().
                                         If None, will try to load from analysis_dir.
        analysis_dir (str, optional): Directory containing analysis results.
        output_dir (str, optional): Directory to save the plots.
                                 If None, will use analysis_dir/plots.

    Returns:
        dict: Dictionary with paths to the generated plots
    """
    # Handle input parameters
    if analysis_results is None and analysis_dir is None:
        logger.error("Either analysis_results or analysis_dir must be provided")
        return {"success": False, "error": "No analysis results provided"}

    # If analysis_results is not provided, try to load from analysis_dir
    if analysis_results is None:
        try:
            # Load summary statistics
            with open(os.path.join(analysis_dir, "summary_statistics.json"), "r") as f:
                summary_statistics = json.load(f)

            # Load correlations
            with open(os.path.join(analysis_dir, "correlations.json"), "r") as f:
                correlations = json.load(f)

            # Load average correlations
            with open(
                os.path.join(analysis_dir, "average_correlations.json"), "r"
            ) as f:
                average_correlations = json.load(f)

            # Load deviations
            with open(os.path.join(analysis_dir, "deviations.json"), "r") as f:
                deviations = json.load(f)

            # Load standardized scores
            standardized_scores = pd.read_csv(
                os.path.join(analysis_dir, "standardized_scores.csv")
            )

            analysis_results = {
                "summary_statistics": summary_statistics,
                "correlations": correlations,
                "average_correlations": average_correlations,
                "deviations": deviations,
                "trait_id": (
                    os.path.basename(os.path.dirname(analysis_dir))
                    if os.path.dirname(analysis_dir)
                    else None
                ),
            }

        except Exception as e:
            logger.error(f"Error loading analysis results: {e}")
            return {"success": False, "error": f"Error loading analysis results: {e}"}
    else:
        # Ensure we have the standardized scores
        if analysis_dir:
            standardized_scores = pd.read_csv(
                os.path.join(analysis_dir, "standardized_scores.csv")
            )
        else:
            logger.error(
                "analysis_dir must be provided if standardized_scores are not in analysis_results"
            )
            return {"success": False, "error": "No standardized scores provided"}

    # Set up output directory
    if output_dir is None:
        if analysis_dir:
            output_dir = os.path.join(analysis_dir, "plots")
        else:
            output_dir = os.path.join(os.getcwd(), "plots")

    os.makedirs(output_dir, exist_ok=True)

    # Get trait name if available
    trait_id = analysis_results.get("trait_id")
    trait_name = None

    if trait_id:
        # Try to get trait name from the trait_id
        try:
            response = requests.get(
                f"https://www.pgscatalog.org/rest/trait/{trait_id}"
            ).json()
            trait_name = response.get("label")
        except:
            pass

    # Generate plots
    plots = {}

    # Ancestry groups
    ancestry_groups = sorted(
        list(
            set(
                group
                for pgs_id in analysis_results["summary_statistics"]
                for group in analysis_results["summary_statistics"][pgs_id]
            )
        )
    )

    # 1. Distribution plots
    for pgs_id in analysis_results["summary_statistics"]:
        plots[f"{pgs_id}_distribution"] = plot_distribution_by_ancestry(
            standardized_scores, pgs_id, output_dir, trait_name
        )

        plots[f"{pgs_id}_standardized"] = plot_standardized_distribution_by_ancestry(
            standardized_scores, pgs_id, output_dir, trait_name
        )

    for group in ancestry_groups:
        plots[f"{group}_distribution"] = plot_distribution_by_pgs(
            standardized_scores, group, output_dir, trait_name
        )

        plots[f"{group}_correlation"] = plot_correlation_matrix(
            standardized_scores, group, output_dir, trait_name
        )

    # 2. Correlation plots
    plots["average_correlations"] = plot_average_correlations(
        analysis_results["average_correlations"], output_dir, trait_name
    )

    # 3. Deviation plots
    plots["deviations"] = plot_deviations(
        analysis_results["deviations"], ancestry_groups, output_dir, trait_name
    )

    return {"success": True, "plots": plots}
