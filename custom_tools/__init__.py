"""
custom_tools: personal bioinformatics/genomics helper library.

Deliberately does not import any submodules here. `muon_preprocessing.py` alone
pulls in muon/scanpy/anndata/mudata/pysam/networkx, and the various subpackages
pull in their own optional heavy dependencies (biopython, pybiomart, mygene,
pybedtools, duckdb). This library is used across many different project
environments that don't all have every optional dependency installed, so
eagerly importing everything here would make a plain `import custom_tools`
fail in any environment missing even one of them.

Import directly from the module/subpackage you need instead, e.g.:

    from custom_tools.peak_utils import format_peak_dataframe
    from custom_tools.download_utils import download_genome_fasta
    from custom_tools.plotting_utils import plot_auroc_auprc
"""

__version__ = "0.0.1"
