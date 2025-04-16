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
import requests

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


def add_value_labels(bars, precision=3, fontsize=None):
    """
    Add value labels on top of each bar in a bar plot.

    Args:
        bars (matplotlib.container.BarContainer): The bar container returned by plt.bar()
        precision (int): Number of decimal places to show
        fontsize (int, optional): Font size for the value labels

    Returns:
        None
    """
    for bar in bars:
        height = bar.get_height()
        if not np.isnan(height):
            plt.text(
                bar.get_x() + bar.get_width() / 2.0,
                height,
                f"{height:.{precision}f}",
                ha="center",
                va="bottom",
                fontsize=fontsize,
            )


def create_bar_plot(
    x_values,
    y_values,
    show_error_bars=False,
    error_values=None,
    bar_width=None,
    color=None,
    alpha=1,
    label=None,
):
    """
    Create a bar plot with optional error bars.

    Args:
        x_values (list): X values for the bars
        y_values (list): Height values for the bars
        show_error_bars (bool): Whether to show error bars
        error_values (list, optional): Error values for error bars
        bar_width (float, optional): Width of the bars
        color: Color of the bars
        alpha (float): Alpha transparency for the bars
        label (str, optional): Label for the bars in the legend

    Returns:
        matplotlib.container.BarContainer: The bar container
    """
    # Create bar plot
    if bar_width is not None:
        bars = plt.bar(
            x_values, y_values, bar_width, alpha=alpha, color=color, label=label
        )
    else:
        bars = plt.bar(x_values, y_values, alpha=alpha, color=color, label=label)

    # Add value labels
    add_value_labels(bars)

    # Add error bars if requested
    if show_error_bars and error_values is not None:
        plt.errorbar(
            x_values,
            y_values,
            yerr=error_values,
            fmt="none",
            capsize=5,
            color="black",
        )

    return bars


def save_plot(output_dir, subdir, filename, close_plot=True):
    """
    Save plot to a file.

    Args:
        output_dir (str): Base output directory
        subdir (str): Subdirectory within the output directory
        filename (str): Filename for the plot
        close_plot (bool): Whether to close the plot after saving

    Returns:
        str: Path to the saved plot, or None if not saved
    """
    if output_dir:
        # Create directory if it doesn't exist
        full_dir = os.path.join(output_dir, subdir)
        os.makedirs(full_dir, exist_ok=True)

        # Save plot
        plot_path = os.path.join(full_dir, filename)
        plt.savefig(plot_path)

        if close_plot:
            plt.close()

        return plot_path
    else:
        plt.show()

        if close_plot:
            plt.close()

        return None


def plot_kde(scores, label, alpha=0.8):
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

    # Get unique ancestry groups and sort alphabetically
    ancestry_groups = sorted(list(scores_df["GROUP"].unique()) + ["ALL"])

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
    return save_plot(output_dir, "distributions", f"{pgs_id}_distributions.png")


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
    return save_plot(output_dir, "distributions", f"{group}_distributions.png")


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

    # Get unique ancestry groups and sort alphabetically
    ancestry_groups = sorted(list(scores_df["GROUP"].unique()) + ["ALL"])

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
        plot_kde(group_z_scores, group)

    # Set up plot
    setup_plot(
        f"Standardized PGS Scores for {pgs_id}",
        "Z-Score",
        "Density",
        add_zero_line=True,
        trait_name=trait_name,
    )

    # Save plot
    return save_plot(
        output_dir,
        "standardized_distributions",
        f"{pgs_id}_standardized_distributions.png",
    )


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
        return save_plot(output_dir, "correlations", f"{group}_correlation_matrix.png")

    except Exception as e:
        logger.warning(f"Error calculating correlation matrix for group {group}: {e}")
        plt.close()
        return None


def plot_average_correlations(
    average_correlations, output_dir=None, trait_name=None, show_error_bars=False
):
    """
    Plot average correlations across ancestry groups with optional error bars.

    Args:
        average_correlations (dict): Dictionary with average correlations by ancestry group
        output_dir (str, optional): Directory to save the plot
        trait_name (str, optional): Trait name to add to title
        show_error_bars (bool): Whether to show error bars. Default is False.

    Returns:
        str or None: Path to the saved plot, or None if no plot was created
    """
    plt.figure(figsize=(12, 7))

    # Sort ancestry groups alphabetically
    sorted_groups = sorted(average_correlations.keys())

    # Extract mean and std for each group
    means = [average_correlations[group]["mean"] for group in sorted_groups]
    stds = [average_correlations[group]["std"] for group in sorted_groups]

    # Create bar plot with optional error bars
    create_bar_plot(sorted_groups, means, show_error_bars, stds)

    # Set up plot
    plt.xlabel("Ancestry Group")
    plt.ylabel("Average Correlation")

    title = "Average Correlation of PGS Scores by Ancestry Group"
    if trait_name:
        title += f" ({trait_name})"

    plt.title(title)
    plt.ylim(bottom=0)
    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()

    # Save plot
    return save_plot(output_dir, "correlations", "average_correlation_bar_chart.png")


def plot_variance_by_ancestry(
    variance_results, output_dir=None, trait_name=None, show_error_bars=False
):
    """
    Plot average variance of z-scores across PGS studies for each ancestry group.

    Lower variance indicates more stable predictions across different PGS models.

    Args:
        variance_results (dict): Dictionary with variance information by ancestry group
        output_dir (str, optional): Directory to save the plot
        trait_name (str, optional): Trait name to add to title
        show_error_bars (bool): Whether to show error bars. Default is False.

    Returns:
        str or None: Path to the saved plot, or None if no plot was created
    """
    plt.figure(figsize=(12, 7))

    # Sort groups alphabetically
    sorted_groups = sorted(variance_results.keys())

    # Extract average variance for each group
    avg_variances = [
        variance_results[group]["average_variance"] for group in sorted_groups
    ]

    # Extract standard deviation if needed for error bars
    std_variances = (
        [variance_results[group]["std_variance"] for group in sorted_groups]
        if show_error_bars
        else None
    )

    # Create bar plot with optional error bars
    create_bar_plot(sorted_groups, avg_variances, show_error_bars, std_variances)

    # Set up plot
    plt.xlabel("Ancestry Group")
    plt.ylabel("Average Variance of Z-Scores")

    title = "Average Variance of Z-Scores Across PGS Studies by Ancestry Group"
    if trait_name:
        title += f" ({trait_name})"

    plt.title(title)
    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()

    # Save plot
    return save_plot(output_dir, "variance", "average_z_score_variance_by_ancestry.png")


def plot_pgs_variance(
    pgs_variance, output_dir=None, trait_name=None, show_error_bars=False
):
    """
    Plot variance of each PGS from the "true" z-score (average across all PGS) by ancestry group.

    Lower variance indicates the PGS is more consistent with the consensus prediction.

    Args:
        pgs_variance (dict): Dictionary with PGS variance information by ancestry group
        output_dir (str, optional): Directory to save the plot
        trait_name (str, optional): Trait name to add to title
        show_error_bars (bool): Whether to show error bars. Default is False.

    Returns:
        str or None: Path to the saved plot, or None if no plot was created
    """
    plt.figure(figsize=(14, 8))

    # Get all unique PGS IDs across all groups
    all_pgs = set()
    for group in pgs_variance:
        all_pgs.update(pgs_variance[group].keys())
    all_pgs = sorted(list(all_pgs))

    # Get all ancestry groups
    groups = sorted(pgs_variance.keys())

    # Set up the bar chart
    bar_width = 0.8 / len(groups)
    index = np.arange(len(all_pgs))

    # Set up colors
    colors = plt.cm.tab10(np.linspace(0, 1, len(groups)))

    # Plot bars for each group
    for i, group in enumerate(groups):
        # Extract variances for this group
        variances = []
        for pgs in all_pgs:
            if pgs in pgs_variance[group]:
                variances.append(pgs_variance[group][pgs]["variance"])
            else:
                variances.append(np.nan)

        # Plot bars
        bars = plt.bar(
            index + i * bar_width,
            variances,
            bar_width,
            color=colors[i],
            label=group,
        )

    # Set up labels and title
    plt.xlabel("Polygenic Score")
    plt.ylabel("Variance from Consensus Z-Score")

    title = "PGS Variance from Consensus Z-Score by Ancestry Group"
    if trait_name:
        title += f" ({trait_name})"
    plt.title(title)

    plt.xticks(
        index + bar_width * (len(groups) - 1) / 2, all_pgs, rotation=45, ha="right"
    )
    plt.legend()
    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()

    # Save plot
    return save_plot(output_dir, "variance", "pgs_variance_by_ancestry.png")


def plot_individual_pgs_variance(
    pgs_variance, output_dir=None, trait_name=None, show_error_bars=False
):
    """
    Create individual plots for each PGS showing its variance from the consensus
    z-score across different ancestry groups.

    Args:
        pgs_variance (dict): Dictionary with PGS variance information by ancestry group
        output_dir (str, optional): Directory to save the plots
        trait_name (str, optional): Trait name to add to title
        show_error_bars (bool): Whether to show error bars. Default is False.

    Returns:
        dict: Dictionary mapping PGS IDs to paths of saved plots
    """
    # Get all unique PGS IDs across all groups
    all_pgs = set()
    for group in pgs_variance:
        all_pgs.update(pgs_variance[group].keys())
    all_pgs = sorted(list(all_pgs))

    # Get all ancestry groups
    groups = sorted(pgs_variance.keys())

    # Dictionary to store paths of saved plots
    plot_paths = {}

    # Create a plot for each PGS
    for pgs in all_pgs:
        plt.figure(figsize=(10, 6))

        # Extract variances for this PGS across groups
        variances = []
        for group in groups:
            if pgs in pgs_variance[group]:
                variances.append(pgs_variance[group][pgs]["variance"])
            else:
                variances.append(np.nan)

        # Create bar plot
        bars = plt.bar(groups, variances)

        # Add value labels on top of each bar
        add_value_labels(bars)

        # Set up labels and title
        plt.xlabel("Ancestry Group")
        plt.ylabel("Variance from Consensus Z-Score")

        title = f"Variance of {pgs} from Consensus Z-Score by Ancestry Group"
        if trait_name:
            title += f" ({trait_name})"
        plt.title(title)

        plt.grid(axis="y", alpha=0.3)
        plt.tight_layout()

        # Save plot
        if output_dir:
            plot_path = save_plot(
                output_dir,
                os.path.join("variance", "individual_pgs"),
                f"{pgs}_variance.png",
            )
            plot_paths[pgs] = plot_path
        else:
            plt.show()
            plt.close()

    return plot_paths


def plot_average_pgs_variance_by_group(
    pgs_variance, output_dir=None, trait_name=None, show_error_bars=False
):
    """
    Plot the average PGS variance across all PGS for each ancestry group.

    Args:
        pgs_variance (dict): Dictionary with PGS variance information by ancestry group
        output_dir (str, optional): Directory to save the plot
        trait_name (str, optional): Trait name to add to title
        show_error_bars (bool): Whether to show error bars. Default is False.

    Returns:
        str or None: Path to the saved plot, or None if no plot was created
    """
    plt.figure(figsize=(10, 6))

    # Calculate average variance for each group
    groups = sorted(pgs_variance.keys())
    avg_variances = []
    std_variances = []

    for group in groups:
        # Get all variances for this group
        variances = [
            pgs_variance[group][pgs]["variance"] for pgs in pgs_variance[group]
        ]

        # Calculate mean and standard deviation
        avg_variances.append(np.mean(variances))
        std_variances.append(np.std(variances))

    # Create bar plot with optional error bars
    create_bar_plot(groups, avg_variances, show_error_bars, std_variances)

    # Set up labels and title
    plt.xlabel("Ancestry Group")
    plt.ylabel("Average Variance from Consensus Z-Score")

    title = "Average PGS Variance from Consensus Z-Score by Ancestry Group"
    if trait_name:
        title += f" ({trait_name})"
    plt.title(title)

    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()

    # Save plot
    return save_plot(output_dir, "variance", "average_pgs_variance_by_group.png")


def visualize_analysis(
    analysis_results=None, analysis_dir=None, output_dir=None, show_error_bars=False
):
    """
    Visualize PGS analysis results.

    Args:
        analysis_results (dict, optional): Analysis results from analyze_scores().
                                         If None, will try to load from analysis_dir.
        analysis_dir (str, optional): Directory containing analysis results.
        output_dir (str, optional): Directory to save the plots.
                                 If None, will use analysis_dir/plots.
        show_error_bars (bool): Whether to show error bars for all plots. Default is False.

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

            # Load individual variance (previously just "variance")
            with open(os.path.join(analysis_dir, "individual_variance.json"), "r") as f:
                individual_variance = json.load(f)

            # Load PGS variance
            with open(os.path.join(analysis_dir, "pgs_variance.json"), "r") as f:
                pgs_variance = json.load(f)

            # Load standardized scores
            standardized_scores = pd.read_csv(
                os.path.join(analysis_dir, "standardized_scores.csv")
            )

            analysis_results = {
                "summary_statistics": summary_statistics,
                "correlations": correlations,
                "average_correlations": average_correlations,
                "individual_variance": individual_variance,
                "pgs_variance": pgs_variance,
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
            trait_name = response["label"]
            trait_name = " ".join(word.capitalize() for word in trait_name.split())
        except:
            logger.warning(f"Could not find trait name for {trait_id}")
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
        analysis_results["average_correlations"],
        output_dir,
        trait_name,
        show_error_bars=show_error_bars,
    )

    # 3. Variance plot
    plots["variance_by_ancestry"] = plot_variance_by_ancestry(
        analysis_results["individual_variance"],
        output_dir,
        trait_name,
        show_error_bars=show_error_bars,
    )

    # 4. PGS variance plot
    plots["pgs_variance"] = plot_pgs_variance(
        analysis_results["pgs_variance"],
        output_dir,
        trait_name,
        show_error_bars=show_error_bars,
    )

    # 5. Individual PGS variance plots
    individual_plots = plot_individual_pgs_variance(
        analysis_results["pgs_variance"],
        output_dir,
        trait_name,
        show_error_bars=show_error_bars,
    )
    plots["individual_pgs_variance"] = individual_plots

    # 6. Average PGS variance by group plot
    plots["average_pgs_variance_by_group"] = plot_average_pgs_variance_by_group(
        analysis_results["pgs_variance"],
        output_dir,
        trait_name,
        show_error_bars=show_error_bars,
    )

    return {"success": True, "plots": plots}
