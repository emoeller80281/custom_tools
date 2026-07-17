"""
Lazy re-exports for custom_tools.download_utils.

    from custom_tools.download_utils import download_genome_fasta  # or any name below

Unlike plotting_utils/stat_utils/standardization_utils, each download_*.py
file here depends on a *different* optional third-party package
(download_genomic_data -> pybiomart+pysam, download_protein_fasta -> biopython,
download_chip_atlas -> duckdb, ...). Eagerly running `from .x import y` for
all of them at package-import time would mean `import custom_tools.download_utils`
fails outright if even one of those packages is missing from the current
environment — even if you only wanted a function from a different file.

Module `__getattr__` (PEP 562) defers each submodule's import until one of its
names is actually accessed, so `download_string_v12_files` works fine even in
an environment that has never installed biopython.
"""

from importlib import import_module

_EXPORTS = {
    "stream_download": ".common_download_utils",
    "download_gene_tss_file": ".download_genomic_data",
    "download_genome_fasta": ".download_genomic_data",
    "download_chrom_sizes": ".download_genomic_data",
    "download_ncbi_gene_info": ".download_genomic_data",
    "download_ensembl_gtf": ".download_genomic_data",
    "fetch_chip_atlas_tf": ".download_chip_atlas",
    "fetch_chip_atlas_tf_list_to_parquet": ".download_chip_atlas",
    "build_chip_atlas_df_from_parquet": ".download_chip_atlas",
    "create_organism_chip_atlas_file": ".download_chip_atlas",
    "download_jaspar_pfms": ".download_jaspar_pfms",
    "download_gene_protein_fastas": ".download_protein_fasta",
    "download_string_v12_files": ".download_string_pdb",
}

__all__ = list(_EXPORTS)


def __getattr__(name: str):
    if name in _EXPORTS:
        module = import_module(_EXPORTS[name], __name__)
        return getattr(module, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(__all__)
