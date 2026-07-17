"""
Re-exports the gene canonicalizer:

    from custom_tools.standardization_utils import GeneCanonicalizer

Only one submodule exists here, and `GeneCanonicalizer` is the entire point of
importing this subpackage, so there's no partial-failure scenario to guard
against by deferring the import (if `mygene` isn't installed, you need to know
that immediately, not the first time you call a method).
"""

from .gene_canonicalizer import GeneCanonicalizer

__all__ = [
    "GeneCanonicalizer",
]
