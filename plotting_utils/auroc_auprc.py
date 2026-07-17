from typing import Literal
import numpy as np
import matplotlib.pyplot as plt

from sklearn.metrics import (
    roc_auc_score,
    average_precision_score,
    precision_recall_curve, 
    roc_curve,
)

def _create_random_distribution(scores, seed: int = 42) -> np.ndarray:
    """
    Create a random distribution of scores with the same shape as the input scores.
    """
    rng = np.random.default_rng(seed)
    arr = np.asarray(scores)   # works for Series or ndarray, no copy if already ndarray
    return rng.uniform(arr.min(), arr.max(), size=arr.shape[0])


def plot_auroc_auprc(
    labels,
    scores,
    roc_line_color="#4195df",
    prc_line_color="#4195df",
    rand_line_color="#747474",
    title=None,
    plot_type: Literal["both", "roc", "prc"] = "both",
) -> plt.Figure:
    """
    labels: array-like of 0/1 labels
    scores: array-like of predicted probabilities after sigmoid

    plot_type:
        "both" -> plot ROC and PRC
        "roc"  -> plot ROC only
        "prc"  -> plot PRC only
    """

    if plot_type not in {"both", "roc", "prc"}:
        raise ValueError("plot_type must be one of: 'both', 'roc', or 'prc'.")

    labels = np.asarray(labels).astype(int).ravel()
    scores = np.asarray(scores).astype(float).ravel()

    rand_scores = _create_random_distribution(scores)

    plots_to_make = []
    if plot_type in {"both", "roc"}:
        plots_to_make.append("roc")
    if plot_type in {"both", "prc"}:
        plots_to_make.append("prc")

    ncols = len(plots_to_make)
    figsize = (7, 4.8) if ncols == 2 else (4.5, 4.8)

    fig, axes = plt.subplots(
        nrows=1,
        ncols=ncols,
        figsize=figsize,
    )

    if ncols == 1:
        axes = [axes]

    for ax, current_plot in zip(axes, plots_to_make):
        
        ax.set_box_aspect(1)
        
        if current_plot == "roc":
            auroc = roc_auc_score(labels, scores)
            fpr, tpr, _ = roc_curve(labels, scores)

            rand_fpr, rand_tpr, _ = roc_curve(labels, rand_scores)
            rand_auroc = roc_auc_score(labels, rand_scores)

            roc_line, = ax.plot(
                fpr,
                tpr,
                lw=2,
                color=roc_line_color,
                label=f"AUROC = {auroc:.3f}",
                zorder=3,
            )

            rand_roc_line, = ax.plot(
                rand_fpr,
                rand_tpr,
                color=rand_line_color,
                linestyle="--",
                lw=2,
                label=f"Random = {rand_auroc:.3f}",
                zorder=2,
            )

            ax.plot(
                [0, 1],
                [0, 1],
                "k--",
                lw=1,
                alpha=0.5,
                zorder=1,
            )

            ax.set_xlabel("False Positive Rate", fontsize=12)
            ax.set_ylabel("True Positive Rate", fontsize=12)
            ax.set_title("AUROC", fontsize=12)

            ax.legend(
                handles=[roc_line, rand_roc_line],
                labels=[f"AUROC = {auroc:.3f}", f"Random = {rand_auroc:.3f}"],
                bbox_to_anchor=(0.5, -0.13),
                loc="upper center",
                borderaxespad=0.0,
                facecolor="none",
                edgecolor="none",
            )

            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)

        elif current_plot == "prc":
            auprc = average_precision_score(labels, scores)
            prec, rec, _ = precision_recall_curve(labels, scores)

            rand_prec, rand_rec, _ = precision_recall_curve(labels, rand_scores)
            rand_auprc = average_precision_score(labels, rand_scores)

            pr_line, = ax.plot(
                rec,
                prec,
                lw=2,
                color=prc_line_color,
                label=f"AUPRC = {auprc:.3f}",
                zorder=3,
            )

            rand_pr_line, = ax.plot(
                rand_rec,
                rand_prec,
                color=rand_line_color,
                linestyle="--",
                lw=2,
                label=f"Random = {rand_auprc:.3f}",
                zorder=2,
            )

            ax.set_xlabel("Recall", fontsize=12)
            ax.set_ylabel("Precision", fontsize=12)
            ax.set_title("AUPRC", fontsize=12)

            ax.legend(
                handles=[pr_line, rand_pr_line],
                labels=[f"AUPRC = {auprc:.3f}", f"Random = {rand_auprc:.3f}"],
                bbox_to_anchor=(0.5, -0.28),
                loc="upper center",
                borderaxespad=0.0,
            )

            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)

    if title is not None:
        fig.suptitle(title, fontsize=12)

    fig.tight_layout()

    return fig