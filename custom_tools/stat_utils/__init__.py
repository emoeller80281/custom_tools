"""
Re-exports the classification metric helper:

    from custom_tools.stat_utils import compute_binary_classification_metrics

Only one submodule exists here and its only dependencies are numpy/sklearn,
so there's no eager-import risk to guard against.
"""

from .classification_metrics import compute_binary_classification_metrics

__all__ = [
    "compute_binary_classification_metrics",
]
