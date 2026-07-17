# custom_tools

Personal library of reusable bioinformatics/genomics helper functions (Eric Moeller), mostly built
around single-cell multiomics (scRNA + scATAC), TF-peak-gene regulatory network construction, and
supporting genomic data download/formatting. Sphinx docs are auto-built from these docstrings
(`docs/`, published via GitHub Actions) — **write new functions with NumPy-style docstrings** to stay
consistent with the rest of the package.

There is no `pyproject.toml`/`setup.py` at the repo root — it's used by adding the repo root to
`sys.path`/`PYTHONPATH` (see how `docs/source/conf.py` does `sys.path.insert(0, ...)`), not via `pip install`.

### Import style

`plotting_utils`, `stat_utils`, and `standardization_utils` re-export their public functions/classes
in `__init__.py`, so either of these works:

```python
from custom_tools.plotting_utils import plot_auroc_auprc            # via subpackage re-export
from custom_tools.plotting_utils.auroc_auprc import plot_auroc_auprc  # direct submodule import, always works
```

`download_utils/__init__.py` also re-exports its functions, but lazily via module `__getattr__`
(PEP 562) rather than plain `from .x import y` — each `download_*.py` file needs a different optional
third-party package (pybiomart, biopython, pysam, duckdb, ...), so an eager re-export would make
`import custom_tools.download_utils` fail outright if even one of those isn't installed, even when you
only wanted a function from a different file. Accessing e.g. `download_utils.download_genome_fasta`
imports just `download_genomic_data.py` at that point, not the others.

The top-level `custom_tools/__init__.py` intentionally re-exports nothing — `muon_preprocessing.py`
alone pulls in muon/scanpy/anndata/pysam, and this library is used across many project environments
that don't all have every optional dependency installed. `peak_utils.py`, `onehot_dna.py`, and
`muon_preprocessing.py` are plain top-level modules (not subpackages), so they were never wrapped by a
re-export layer to begin with — always import from them directly, e.g. `from custom_tools.peak_utils
import format_peak_dataframe`.

## Where to look for what

| Need | Module |
|---|---|
| Peak ID parsing/formatting (`chr1:100-200` normalization) | `peak_utils.py`, `onehot_dna.py` (has its own private copy) |
| Nearest-gene/TSS distance for peaks | `peak_utils.py` |
| One-hot encoding DNA sequence for ML models | `onehot_dna.py` |
| scRNA/scATAC QC, filtering, PCA/UMAP, MOFA+ integration, metacells | `muon_preprocessing.py` |
| Downloading reference genomes, GTFs, TSS, chrom sizes, gene info | `download_utils/download_genomic_data.py` |
| Downloading ChIP-Atlas peaks (ground truth TF-DNA edges) | `download_utils/download_chip_atlas.py` |
| Downloading JASPAR motif PFMs | `download_utils/download_jaspar_pfms.py` |
| Downloading STRING protein-protein interaction files | `download_utils/download_string_pdb.py` |
| Downloading RefSeq protein FASTAs (NCBI Entrez) | `download_utils/download_protein_fasta.py` |
| Generic streaming file download w/ progress bar | `download_utils/common_download_utils.py` |
| Binary classification metrics (AUROC/AUPRC/F1/precision@k/etc.) | `stat_utils/classification_metrics.py` |
| ROC/PR curve plotting | `plotting_utils/auroc_auprc.py` |
| True/false score distribution histograms/KDEs | `plotting_utils/histograms.py` |
| Canonicalizing gene symbols/aliases/Ensembl/Entrez IDs to official symbol | `standardization_utils/gene_canonicalizer.py` |

## `peak_utils.py`

Utilities for genomic peak string handling and peak-to-gene assignment (used across scATAC-seq /
regulatory network pipelines). Peaks are represented as `chrN:start-end` strings.

- `format_individual_peak(peak_id: str) -> str` — normalize a single peak ID from `chrN-start-end`
  or `chrN:start:end` into the canonical `chrN:start-end` form.
- `format_peak_dataframe(peak_ids: pd.Series | pd.Index) -> pd.DataFrame` — vectorized version;
  parses a whole column of peak IDs (handles dashes, colons, `hg38.` prefixes, BED-like rows) into a
  DataFrame with `chromosome`, `start`, `end`, `strand`, `peak_id` columns. Raises on malformed IDs.
  **Prefer this over writing custom peak-string regexes.**
- `get_peak_length(peak_id_col: pd.Series) -> pd.Series` — base-pair length of each peak.
- `find_genes_near_peaks(peak_bed, tss_bed, tss_distance_cutoff=1e6) -> pd.DataFrame` — uses
  `pybedtools.BedTool.window` to find TSS within a cutoff distance of each peak; returns
  `peak_id`/`gene_id`/`TSS_dist`, sorted by distance.
- `set_tg_as_closest_gene_tss_to_peak(tf_peak_edge_df, peaks_gene_distance_file) -> pd.DataFrame` —
  assigns each TF-peak edge its target gene (TG) as the closest gene by TSS distance, reading a
  precomputed peak→gene distance parquet.

## `onehot_dna.py`

Fast DNA sequence one-hot encoding for ML input pipelines (e.g. sequence-based TF binding models).

- `onehot_dna_sequence(seq: str) -> np.ndarray` — vectorized A/C/G/T (case-insensitive) one-hot
  encode, shape `(L, 4)` float32. Uses a lookup table over raw bytes, much faster than naive loops.
- `parse_peak(peak: str) -> tuple[str, int, int]` — parses `chrom:start-end` into components.
- `load_peak_sequence(genome_fasta, selected_peak) -> str` — pulls the sequence for one peak from an
  indexed FASTA via `pyfaidx`.
- `load_chrom_sizes(chromsizes_file) -> dict[str, int]` — parses a UCSC `.chrom.sizes` file.
- `create_centered_peak_onehot_array(peak_ids, genome_fasta, chrom_sizes, peak_id_to_idx, flank_size, ...)`
  — the main entry point: builds a stacked `(n_peaks, 2*flank_size, 4)` one-hot array of sequences
  centered on each peak midpoint, with out-of-bounds N-padding. Supports `num_workers > 1` for
  multiprocess encoding (chunks peaks across a `ProcessPoolExecutor`) and a `tqdm` progress bar.
  Use this rather than hand-rolling a peak→sequence→one-hot loop.

## `muon_preprocessing.py`

End-to-end single-cell multiome (RNA + ATAC) preprocessing pipeline built on `muon`/`scanpy`/`anndata`.
Also runnable as a CLI script (`if __name__ == "__main__"`, driven by `argparse` + a
`qc_filtering_settings.tsv` per-sample threshold table). Best used as a reference/template for a new
multiome QC pipeline rather than imported piecemeal, but individual functions are still useful standalone:

- `filter_to_human(mdata)` — for barnyard (mixed-species) experiments, keep only `hg38`-prefixed ATAC
  peaks and strip the prefix.
- `create_fragment_index_file(frag_path)` — tabix-index a `fragments.tsv.gz` (idempotent).
- `construct_mdata_from_gene_by_cell_matrices(rna_count_file, atac_count_file) -> mu.MuData` — build a
  MuData object directly from gene×cell / peak×cell CSVs (alternative to 10x mtx/h5 input).
- `load_raw_data(sample_name, sample_data_dir, ...)` — auto-detects and loads 10x mtx triplets, a raw
  `.h5`/`.h5mu` file, or raw count CSVs for a sample directory.
- `normalize_peak_format(peak_id)` — same peak-string normalization as `peak_utils.format_individual_peak`
  (duplicated here; prefer `peak_utils` version for new code).
- `MudataProcessor` class — stateful QC/preprocessing pipeline object (`self.mdata`, `self.rna`,
  `self.atac`):
  - `rna_qc_filter(...)` — mito filtering, gene/count thresholds, TF-gene allowlist retention,
    normalization, log1p, HVG selection, scaling; writes RNA QC violin plots if `fig_dir` given.
  - `rna_pca_and_neighbors(...)` — PCA, neighbors, UMAP, Leiden clustering for RNA.
  - `construct_peak_annotation(...)` — assigns each ATAC peak a `promoter`/`distal`/`intergenic` label
    and nearest gene based on TSS proximity (10x-style annotation TSV).
  - `atac_qc_filter(...)` — peak/count thresholds, TF-IDF, LSI, HVG peaks, neighbors/PCA/UMAP/Leiden,
    calls `construct_peak_annotation` internally.
  - `nucleosome_signal(frag_path, ...)`, `tss_enrichment(frag_path, ...)`,
    `tss_enrichment_plot(...)` — ATAC fragment-based QC metrics/plots.
  - `save_stability_subsamplings(...)`, `save_mdata()`, `save_rna()`, `save_atac()` — persistence
    helpers.
- `integrate_rna_atac(mdata, sample_processed_data_dir, sample_name, fig_dir=None)` — joint RNA+ATAC
  integration via MOFA+ (`mu.tl.mofa`), with guards against non-finite/zero-variance features that
  would otherwise desync MOFA's feature dimensions; produces joint UMAP/Leiden clusters and
  rank_genes/rank_peaks group markers.
- `save_processed_data(mdata, sample_processed_data_dir)` — export RNA/ATAC as feature×cell parquet
  (`scRNA_seq_processed.parquet` / `scATAC_seq_processed.parquet`) plus the full `.h5mu`.
- `create_metacells(mdata, sample_processed_data_dir, hops=2)` — builds neighbor-graph diffusion
  "metacell" pseudobulk profiles (row-normalized adjacency matrix powered by `hops`) and saves
  `TG_pseudobulk.parquet` / `RE_pseudobulk.parquet`. Useful for denoising sparse single-cell data
  before downstream network inference.

## `download_utils/`

- `common_download_utils.stream_download(url, dest, chunk=1<<20, desc=None)` — the shared streaming
  download-with-progress-bar-and-atomic-rename helper; other download modules either call this or
  duplicate its pattern inline. **Prefer calling this directly for any new download code** instead of
  re-writing the requests/tqdm/tmp-file-rename boilerplate.
- `download_genomic_data.py` — reference genome/annotation fetchers for mouse (`mm10`) and human (`hg38`):
  - `download_gene_tss_file(save_file, gene_dataset_name="hsapiens_gene_ensembl", ensembl_version=...)`
    — pulls TSS coordinates from Ensembl BioMart via `pybiomart`, writes a BED-like TSV.
  - `download_genome_fasta(organism_code, save_dir)` — downloads UCSC genome FASTA, converts to
    BGZF (using `bgzip` CLI if available, else pure-Python fallback), and indexes with `pysam.faidx`.
  - `download_chrom_sizes(organism_code, save_dir)` — UCSC `.chrom.sizes`.
  - `download_ncbi_gene_info(organism_code, out_path=None)` — NCBI `gene_info.gz` (used by
    `GeneCanonicalizer.load_ncbi_gene_info`).
  - `download_ensembl_gtf(organism_code, release=None, assembly=None, out_dir=None, decompress=False)`
    — Ensembl GTF (used by `GeneCanonicalizer.load_gtf`).
- `download_chip_atlas.py` — fetch ChIP-Atlas TF ChIP-seq peaks as ground-truth TF-DNA edges:
  `fetch_chip_atlas_tf` (single TF), `fetch_chip_atlas_tf_list_to_parquet` (parallel, ThreadPoolExecutor,
  skips already-downloaded TFs), `build_chip_atlas_df_from_parquet` (DuckDB merge into one parquet),
  `create_organism_chip_atlas_file` (end-to-end: fetch all TFs for a species → combined parquet, cached).
- `download_jaspar_pfms.py` — `download_jaspar_pfms(save_dir, tax_id="10090", version=2024, max_workers=8)`
  pulls all JASPAR motif PFMs for an organism via the REST API (parallel download, skips existing files).
- `download_protein_fasta.py` — `download_gene_protein_fastas(gene_names, organism, output_dir, email, ...)`
  fetches one representative RefSeq protein FASTA per gene via NCBI Entrez (prefers RefSeq Select / NP_
  over XP_/low-quality records), rate-limited with retry/backoff.
- `download_string_pdb.py` — `download_string_v12_files(string_dir, string_org_code)` downloads STRING
  v12.0 protein info + detailed protein links (PPI network) for a taxonomy code.

## `stat_utils/classification_metrics.py`

- `compute_binary_classification_metrics(labels, scores, score_threshold=0.5, random_state=42) -> dict`
  — one-call summary for evaluating any binary classifier/edge-prediction task (e.g. TF-target edge
  prediction): AUROC, AUPRC, accuracy, precision, recall, F1, **early_precision** (precision restricted
  to top 10,000 highest-scored predictions — useful for large sparse network evaluation), plus random-baseline
  AUROC/AUPRC for comparison and `n_pos`/`n_neg`/`n_edges` counts. Returns NaNs gracefully if only one
  class present. **Use this instead of manually calling sklearn metrics one at a time.**

## `plotting_utils/`

- `auroc_auprc.plot_auroc_auprc(labels, scores, plot_type="both"|"roc"|"prc", title=None, ...) -> plt.Figure`
  — side-by-side ROC and/or PR curve plot with a random-baseline curve overlaid for comparison.
- `histograms.plot_score_histograms(labels, scores, panel_kind="kde"|"hist", balance_pos_neg=True, ...) -> plt.Figure`
  — overlaid true/false score distribution (KDE or histogram), auto-downsampling the majority class
  (`_balance_pos_neg`) so the comparison isn't dominated by class imbalance (typical for regulatory
  network edge scores where negatives vastly outnumber positives).

Both plotting functions expect binary `labels` (0/1) and continuous `scores`; pair naturally with
`stat_utils.classification_metrics.compute_binary_classification_metrics` when evaluating a model/edge
scoring method.

## `standardization_utils/gene_canonicalizer.py`

- `GeneCanonicalizer(species="10090", use_mygene=True)` — resolves messy gene identifiers (Excel-mangled
  symbols like `1-Sep`, Greek letters, aliases, Ensembl IDs, Entrez IDs) to one official gene symbol.
  Species-agnostic; taxonomy id string (`"10090"` mouse, `"9606"` human) controls NCBI/MyGene lookups.
  - `load_gtf(gtf_path)` — build Ensembl-ID→symbol map from a GTF (pair with `download_ensembl_gtf`).
  - `load_ncbi_gene_info(gene_info_path, species_taxid="10090")` — build Entrez-ID→symbol and
    alias→symbol maps (pair with `download_ncbi_gene_info`).
  - `canonical_symbol(s) -> str` — resolve one identifier (curated overrides → local GTF/NCBI maps →
    MyGene.info fallback, cached).
  - `canonicalize_series(s: pd.Series, batch_size=5000) -> pd.Series` — vectorized/batched version for
    a whole column; batches unresolved IDs through MyGene.info rather than one-by-one.
  - `standardize_df(df, tf_col, tg_col) -> pd.DataFrame` — canonicalize TF/TG columns of an edge list
    in one call and drop rows that fail to resolve.
  - `coverage_report()` — counts of loaded official symbols / Ensembl / Entrez / alias mappings, useful
    for sanity-checking before trusting canonicalization on a new species/dataset.
  - Typical setup: `download_ensembl_gtf` + `download_ncbi_gene_info` → `load_gtf` + `load_ncbi_gene_info`
    → `standardize_df` on any TF-target edge table before merging datasets from different sources
    (this is exactly the kind of alias mismatch that silently breaks joins across ChIP-Atlas/JASPAR/
    STRING/RNA data sources).

Do not alter the section below when updating this file:

When updating the repository, ensure that you:
 - Update`CLAUDE.md` to keep the documentation up-to-date.
 - Update `docs/requirements.txt` if any packages have been added.
 - Ensure that any new functions have full docstrings and that arguments are typed.