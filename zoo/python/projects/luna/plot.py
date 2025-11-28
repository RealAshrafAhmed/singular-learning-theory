import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde, chi2
from scipy.spatial import ConvexHull
from matplotlib.patches import Ellipse
import seaborn as sns


def visualize_joint_confidence_regions(
    mle_samples,
    true_params=None,
    confidence_levels=[0.68, 0.95],
    methods=["scatter", "kde", "ellipse", "hull"],
    figsize=(15, 12),
):
    """
    Visualize joint confidence regions from Monte Carlo MLE samples.

    Parameters:
    -----------
    mle_samples : list of dicts
        Each dict contains MLE estimates: {'p1': val, 'p2': val, 'rho1': val, 'rho2': val}
    true_params : dict
        True parameter values for reference
    confidence_levels : list
        Confidence levels to plot (e.g., [0.68, 0.95])
    methods : list
        Which methods to show: 'scatter', 'kde', 'ellipse', 'hull'
    """

    # Extract parameter pairs
    param_names = mle_samples.keys()
    n_params = len(param_names)

    # Create parameter arrays
    # params = {}
    # for name in param_names:
    #     params[name] = np.array([sample[name] for sample in mle_samples])

    # Create subplot grid for pairwise plots
    n_pairs = n_params * (n_params - 1) // 2  # Number of unique pairs
    n_methods = len(methods)

    fig, axes = plt.subplots(n_pairs, n_methods, figsize=figsize)
    if n_pairs == 1:
        axes = axes.reshape(1, -1)
    if n_methods == 1:
        axes = axes.reshape(-1, 1)

    colors = ["blue", "red", "green", "orange"]

    pair_idx = 0
    for i in range(n_params):
        for j in range(i + 1, n_params):
            param_x, param_y = param_names[i], param_names[j]
            x_data, y_data = mle_samples[param_x], mle_samples[param_y]

            for method_idx, method in enumerate(methods):
                ax = axes[pair_idx, method_idx]

                if method == "scatter":
                    plot_scatter_method(
                        ax,
                        x_data,
                        y_data,
                        param_x,
                        param_y,
                        true_params,
                        confidence_levels,
                    )
                elif method == "kde":
                    plot_kde_method(
                        ax,
                        x_data,
                        y_data,
                        param_x,
                        param_y,
                        true_params,
                        confidence_levels,
                    )
                elif method == "ellipse":
                    plot_ellipse_method(
                        ax,
                        x_data,
                        y_data,
                        param_x,
                        param_y,
                        true_params,
                        confidence_levels,
                    )
                elif method == "hull":
                    plot_hull_method(
                        ax,
                        x_data,
                        y_data,
                        param_x,
                        param_y,
                        true_params,
                        confidence_levels,
                    )

                if pair_idx == 0:  # Add method title to top row
                    ax.set_title(
                        f"{method.capitalize()} Method", fontsize=12, fontweight="bold"
                    )

                if method_idx == 0:  # Add parameter pair label to leftmost column
                    ax.set_ylabel(
                        f"{param_y} vs {param_x}", fontsize=10, rotation=90, labelpad=20
                    )

            pair_idx += 1

    plt.tight_layout()
    return fig, axes


def plot_scatter_method(
    ax, x_data, y_data, param_x, param_y, true_params, confidence_levels
):
    """Scatter plot with truth point."""
    ax.scatter(
        x_data,
        y_data,
        alpha=0.6,
        s=15,
        color="steelblue",
        label=f"MLE samples (n={len(x_data)})",
    )

    if true_params:
        ax.scatter(
            true_params[param_x],
            true_params[param_y],
            color="red",
            s=100,
            marker="*",
            label="True parameters",
            zorder=10,
        )

    ax.set_xlabel(param_x)
    ax.set_ylabel(param_y)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)


def plot_kde_method(
    ax, x_data, y_data, param_x, param_y, true_params, confidence_levels
):
    """KDE contour plot."""
    ax.scatter(x_data, y_data, alpha=0.4, s=8, color="lightblue")

    try:
        # Create KDE
        kde = gaussian_kde([x_data, y_data])

        # Create grid
        x_min, x_max = x_data.min(), x_data.max()
        y_min, y_max = y_data.min(), y_data.max()

        x_pad = (x_max - x_min) * 0.1
        y_pad = (y_max - y_min) * 0.1

        xx, yy = np.mgrid[
            x_min - x_pad : x_max + x_pad : 50j, y_min - y_pad : y_max + y_pad : 50j
        ]

        positions = np.vstack([xx.ravel(), yy.ravel()])
        density = kde(positions).reshape(xx.shape)

        # Plot contours for each confidence level
        colors = ["orange", "red"]
        for i, conf_level in enumerate(confidence_levels):
            # Find contour level
            density_flat = density.flatten()
            density_sorted = np.sort(density_flat)[::-1]
            threshold_idx = int(conf_level * len(density_sorted))
            threshold = density_sorted[threshold_idx]

            ax.contour(
                xx,
                yy,
                density,
                levels=[threshold],
                colors=[colors[i]],
                linewidths=2,
                alpha=0.8,
            )
            ax.contourf(
                xx,
                yy,
                density,
                levels=[threshold, density.max()],
                colors=[colors[i]],
                alpha=0.2,
            )

    except Exception as e:
        print(f"KDE failed for {param_x} vs {param_y}: {e}")
        ax.text(0.5, 0.5, "KDE Failed", transform=ax.transAxes, ha="center")

    if true_params:
        ax.scatter(
            true_params[param_x],
            true_params[param_y],
            color="red",
            s=100,
            marker="*",
            zorder=10,
        )

    ax.set_xlabel(param_x)
    ax.set_ylabel(param_y)
    ax.grid(True, alpha=0.3)


def plot_ellipse_method(
    ax, x_data, y_data, param_x, param_y, true_params, confidence_levels
):
    """Confidence ellipses assuming Gaussian distribution."""
    ax.scatter(x_data, y_data, alpha=0.4, s=8, color="lightblue")

    try:
        # Compute sample statistics
        mean_x, mean_y = np.mean(x_data), np.mean(y_data)
        cov_matrix = np.cov(x_data, y_data)

        colors = ["orange", "red"]
        for i, conf_level in enumerate(confidence_levels):
            # Chi-square critical value for 2 DOF
            chi2_val = chi2.ppf(conf_level, df=2)

            # Eigenvalues and eigenvectors
            eigenvals, eigenvecs = np.linalg.eigh(cov_matrix)

            # Sort eigenvalues
            idx = eigenvals.argsort()[::-1]
            eigenvals = eigenvals[idx]
            eigenvecs = eigenvecs[:, idx]

            # Ellipse parameters
            width = 2 * np.sqrt(chi2_val * eigenvals[0])
            height = 2 * np.sqrt(chi2_val * eigenvals[1])
            angle = np.degrees(np.arctan2(eigenvecs[1, 0], eigenvecs[0, 0]))

            ellipse = Ellipse(
                (mean_x, mean_y),
                width,
                height,
                angle=angle,
                facecolor=colors[i],
                alpha=0.2,
                edgecolor=colors[i],
                linewidth=2,
                label=f"{int(conf_level * 100)}% CI",
            )
            ax.add_patch(ellipse)

        # Mark sample mean
        ax.scatter(
            mean_x,
            mean_y,
            color="blue",
            s=50,
            marker="x",
            label="Sample mean",
            zorder=10,
        )

    except Exception as e:
        print(f"Ellipse failed for {param_x} vs {param_y}: {e}")
        ax.text(0.5, 0.5, "Ellipse Failed", transform=ax.transAxes, ha="center")

    if true_params:
        ax.scatter(
            true_params[param_x],
            true_params[param_y],
            color="red",
            s=100,
            marker="*",
            zorder=10,
        )

    ax.set_xlabel(param_x)
    ax.set_ylabel(param_y)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)


def plot_hull_method(
    ax, x_data, y_data, param_x, param_y, true_params, confidence_levels
):
    """Convex hull of inner percentiles."""
    ax.scatter(x_data, y_data, alpha=0.4, s=8, color="lightblue")

    try:
        points = np.column_stack([x_data, y_data])
        centroid = np.mean(points, axis=0)
        distances = np.linalg.norm(points - centroid, axis=1)

        colors = ["orange", "red"]
        for i, conf_level in enumerate(confidence_levels):
            # Take inner conf_level percentage of points
            threshold_idx = int(conf_level * len(distances))
            threshold_distance = np.partition(distances, threshold_idx)[threshold_idx]

            inner_points = points[distances <= threshold_distance]

            if len(inner_points) >= 3:  # Need at least 3 points for hull
                hull = ConvexHull(inner_points)

                # Plot hull boundary
                for simplex in hull.simplices:
                    ax.plot(
                        inner_points[simplex, 0],
                        inner_points[simplex, 1],
                        color=colors[i],
                        linewidth=2,
                        alpha=0.8,
                    )

                # Fill hull
                hull_points = inner_points[hull.vertices]
                ax.fill(
                    hull_points[:, 0],
                    hull_points[:, 1],
                    color=colors[i],
                    alpha=0.2,
                    label=f"{int(conf_level * 100)}% Hull",
                )

    except Exception as e:
        print(f"Hull failed for {param_x} vs {param_y}: {e}")
        ax.text(0.5, 0.5, "Hull Failed", transform=ax.transAxes, ha="center")

    if true_params:
        ax.scatter(
            true_params[param_x],
            true_params[param_y],
            color="red",
            s=100,
            marker="*",
            zorder=10,
        )

    ax.set_xlabel(param_x)
    ax.set_ylabel(param_y)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)


# Example usage function
def example_monte_carlo_confidence_regions():
    """Example of how to use the visualization."""

    # Simulate MLE samples (replace with your actual Monte Carlo results)
    np.random.seed(42)
    n_samples = 1000

    # Simulate some correlated MLE estimates around true values
    true_params = {"p1": 0.3, "p2": 0.7, "rho1": 0.4, "rho2": 0.6}

    mle_samples = []
    for _ in range(n_samples):
        # Add some correlation and noise around true values
        noise = np.random.multivariate_normal(
            [0, 0, 0, 0],
            [
                [0.01, 0.005, 0.002, 0],
                [0.005, 0.01, 0, 0.002],
                [0.002, 0, 0.008, 0.004],
                [0, 0.002, 0.004, 0.008],
            ],
        )

        sample = {
            "p1": true_params["p1"] + noise[0],
            "p2": true_params["p2"] + noise[1],
            "rho1": true_params["rho1"] + noise[2],
            "rho2": true_params["rho2"] + noise[3],
        }

        # Ensure valid parameter ranges
        sample["p1"] = np.clip(sample["p1"], 0.01, 0.99)
        sample["p2"] = np.clip(sample["p2"], 0.01, 0.99)
        sample["rho1"] = np.clip(sample["rho1"], 0.01, 0.99)
        sample["rho2"] = np.clip(sample["rho2"], 0.01, 0.99)

        mle_samples.append(sample)

    # Create visualization
    fig, axes = visualize_joint_confidence_regions(
        mle_samples,
        true_params=true_params,
        confidence_levels=[0.68, 0.95],
        methods=["scatter", "kde", "ellipse", "hull"],
    )

    plt.suptitle(
        "Joint Confidence Regions: Monte Carlo MLE Distribution", fontsize=16, y=0.98
    )
    plt.show()

    return fig, axes


def plot_pairwise_scatter(
    p1_samples, p2_samples, rho1_samples, rho2_samples, figsize=(10, 10)
):
    """
    Create pairwise scatter plots using seaborn.

    Parameters:
    -----------
    p1_samples, p2_samples, rho1_samples, rho2_samples : numpy arrays
        Bootstrap samples for each parameter
    figsize : tuple
        Figure size
    """

    # Create DataFrame
    data = pd.DataFrame(
        {"ρ₁": rho1_samples, "ρ₂": rho2_samples, "p₁": p1_samples, "p₂": p2_samples}
    )

    # Create pairplot
    g = sns.pairplot(
        data,
        diag_kind="hist",  # Histograms on diagonal
        plot_kws={"alpha": 0.6, "s": 20},  # Scatter plot settings
        diag_kws={"alpha": 0.7},
    )  # Histogram settings

    # Customize the plot
    g.fig.suptitle("Bootstrap Samples: Pairwise Scatter Plots", fontsize=16, y=1.02)

    return g
