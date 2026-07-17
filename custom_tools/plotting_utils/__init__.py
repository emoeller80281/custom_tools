"""
Re-exports the plotting helpers so callers don't need to know which file a
function lives in:

    from custom_tools.plotting_utils import plot_auroc_auprc, plot_score_histograms

Both submodules only depend on numpy/matplotlib/seaborn/sklearn, which are
always installed together in any environment that uses this subpackage at
all, so eager imports here carry no real risk of an ImportError on an
unrelated, unused function.
"""

from .auroc_auprc import plot_auroc_auprc
from .histograms import plot_score_histograms

__all__ = [
    "plot_auroc_auprc",
    "plot_score_histograms",
]
