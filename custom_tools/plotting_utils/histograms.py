import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def _balance_pos_neg(labels, scores) -> tuple[np.ndarray, np.ndarray]:
    """
    Balance the number of positive and negative examples by downsampling the majority class.
    
    Parameters
    ----------
    labels : array-like
        True binary labels (0 or 1).
    scores : array-like
        Predicted scores or probabilities.
        
    Returns
    -------
    balanced_labels : np.ndarray
        Balanced binary labels.
    balanced_scores : np.ndarray
        Balanced predicted scores.
    """
    true_scores = scores[labels == 1]
    false_scores = scores[labels == 0]
    
    n_pos = len(true_scores)
    n_neg = len(false_scores)
    
    if n_pos == 0 or n_neg == 0:
        raise ValueError("Both positive and negative examples are required for balancing.")
    
    if n_pos > n_neg:
        true_scores_balanced = np.random.choice(true_scores, size=n_neg, replace=False)
        false_scores_balanced = false_scores    
    elif n_neg > n_pos:
        false_scores_balanced = np.random.choice(false_scores, size=n_pos, replace=False)
        true_scores_balanced = true_scores
        
    balanced_labels = np.concatenate([np.ones_like(true_scores_balanced), np.zeros_like(false_scores_balanced)])
    balanced_scores = np.concatenate([true_scores_balanced, false_scores_balanced])
    
    # Shuffle the balanced dataset
    indices = np.arange(len(balanced_labels))
    np.random.shuffle(indices)
    balanced_labels = balanced_labels[indices]
    balanced_scores = balanced_scores[indices]
    
    return balanced_labels, balanced_scores

def plot_score_histograms(
    labels,
    scores, 
    n_bins=75, 
    random_state=42, 
    y_log=False,
    panel_kind="kde",
    density=False,
    title=None,
    y_lim=None,
    x_lim=None,
    balance_pos_neg: bool = True,
):
    """
    Plot histograms or KDEs of predicted scores for true and false labels.
    
    Parameters
    ----------
    labels : array-like
        True binary labels (0 or 1).
    scores : array-like
        Predicted scores or probabilities.
    n_bins : int
        Number of bins for the histogram (default: 75).
    random_state : int
        Random seed for reproducibility (default: 42).
    y_log : bool
        Whether to use a logarithmic scale for the y-axis (default: False).
    panel_kind : str
        Type of plot to generate: "hist" for histogram, "kde" for kernel density estimate (default: "kde").
    density : bool
        Whether to normalize the histogram to form a probability density (default: False).
    title : str, optional
        Title for the figure. If None, no title is set.
    y_lim : tuple, optional
        Limits for the y-axis (min, max). If None, automatic limits are used.
    x_lim : tuple, optional
        Limits for the x-axis (min, max). If None, automatic limits are used.
    balance_pos_neg : bool
        Whether to balance the number of positive and negative examples by downsampling the majority class (default: True).
        
    Returns
    -------
    matplotlib.figure.Figure
        The generated figure object.
    """
    
    
    fig, ax = plt.subplots(
        nrows=1, 
        ncols=1, 
        figsize=(4, 3),
        squeeze=False,
    )

    y = np.asarray(labels).astype(int).ravel()
    s = np.asarray(scores).astype(float).ravel()

    if balance_pos_neg:
        balanced_labels, balanced_scores = _balance_pos_neg(y, s)

        true_vals = balanced_scores[balanced_labels == 1]
        false_vals = balanced_scores[balanced_labels == 0]
    else:
        true_vals = s[y == 1]
        false_vals = s[y == 0]

    min_len = min(len(true_vals), len(false_vals))
    if min_len == 0:
        raise ValueError("Not enough positives/negatives to plot histograms.")

    rng = np.random.default_rng(random_state)
    true_vals = rng.choice(true_vals, size=min_len, replace=False)
    false_vals = rng.choice(false_vals, size=min_len, replace=False)

    combined = np.concatenate([true_vals, false_vals])
    bins = np.linspace(combined.min(), combined.max(), n_bins)

    plot_ax = ax[0, 0]

    if panel_kind == "hist":
        plot_ax.hist(
            false_vals,
            bins=bins,
            alpha=0.6,
            label="False",
            density=density,
        )
        plot_ax.hist(
            true_vals,
            bins=bins,
            alpha=0.6,
            label="True",
            density=density,
        )

        plot_ax.set_title("True vs False Scores", fontsize=12)
        plot_ax.set_xlabel("Score", fontsize=12)
        plot_ax.set_ylabel("Density" if density else "Count", fontsize=12)
        plot_ax.legend(fontsize=9)

    elif panel_kind == "kde":
        sns.kdeplot(
            false_vals,
            ax=plot_ax,
            label="False",
            fill=True,
            common_norm=False,
            bw_adjust=1.0,
            color="#747474",
        )
        sns.kdeplot(
            true_vals,
            ax=plot_ax,
            label="True",
            fill=True,
            common_norm=False,
            bw_adjust=1.0,
            color="#4195df"
        )

        plot_ax.set_title("True vs False Score Density", fontsize=12)
        plot_ax.set_xlabel("Score", fontsize=12)
        plot_ax.set_ylabel("Density", fontsize=12)
        plot_ax.legend(fontsize=9)

    if y_log:
        plot_ax.set_yscale("log")
        plot_ax.set_ylim(bottom=0.1)
        
    if y_lim:
        plot_ax.set_ylim(y_lim[0], y_lim[1])
    if x_lim:
        plot_ax.set_xlim(x_lim[0], x_lim[1])
        
    if title is not None:
        plt.suptitle(title, fontsize=12)

    fig.tight_layout(rect=[0, 0, 1, 0.98])

    return fig