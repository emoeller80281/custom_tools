import pandas as pd
import numpy as np
import pybedtools

def format_individual_peak(peak_id: str) -> str:
    """
    Normalize peak format from chrN-start-end or chrN:start:end to chrN:start-end.
    Handles both formats as input and always outputs chrN:start-end.
    """
    if not isinstance(peak_id, str):
        return peak_id
    
    # Try to parse chr-start-end format (with dashes)
    parts = peak_id.split('-')
    if len(parts) >= 3:
        # Assume format is chr-start-end where chr might have dashes
        # Work backwards: the last two parts are start and end
        try:
            end = int(parts[-1])
            start = int(parts[-2])
            chrom = '-'.join(parts[:-2])  # Everything before the last two parts
            return f"{chrom}:{start}-{end}"
        except (ValueError, IndexError):
            pass
    
    # Already in chr:start-end format or some other format, return as-is
    return peak_id


def format_peak_dataframe(peak_ids: pd.Series|pd.Index) -> pd.DataFrame:
    """
    Splits peaks from `chrN:start-end` or `chrN-start-end` format into a DataFrame.
    
    Creates a dataframe with the following columns:
    1) "peak_id": peakN+1 where N is the index position of the peak
    2) "chromosome": chrN
    3) "start"
    4) "end"
    5) "strand": List of "." values, we dont have strand information for our peaks.
    
    Args:
        peak_ids (pd.Series):
            Series containing the peak locations in "chrN:start-end" or "chrN-start-end" format.
            
    Returns:
        peak_df (pd.DataFrame):
            DataFrame of peak locations in the correct format for Homer and the sliding window method
    """
    if peak_ids.empty:
        raise ValueError("Input peak ID list is empty.")
    
    peak_ids = pd.Series(peak_ids).drop_duplicates().astype(str).str.strip()

    # Keep only canonical chromosomes/records that contain "chr".
    peak_ids = peak_ids[peak_ids.str.contains("chr", regex=False)]
    if peak_ids.empty:
        raise ValueError("No peak IDs containing 'chr' were found after filtering.")

    # Primary parse: supports strings containing a canonical peak token anywhere,
    # such as "chr1:100-200", "hg38.chr1:100-200", or "hg38.chr1-100-200".
    parsed = peak_ids.str.extract(
        r'.*?(?P<chromosome>chr[^\s:-]+)(?::|-)(?P<start>\d+)-(?P<end>\d+)\s*$'
    )

    # Fallback parse: BED-like row strings, e.g. "hg38.chr1 100 200 ...".
    missing = parsed["chromosome"].isnull() | parsed["start"].isnull() | parsed["end"].isnull()
    if missing.any():
        bed_like = peak_ids[missing].str.extract(
            r'^\s*(?P<chromosome>\S+)\s+(?P<start>\d+)\s+(?P<end>\d+)(?:\s|$)'
        )
        parsed.loc[missing, ["chromosome", "start", "end"]] = bed_like[["chromosome", "start", "end"]].values

    # Normalize chromosomes like "hg38.chr1" to "chr1" by removing any prefix before "chr".
    parsed["chromosome"] = parsed["chromosome"].str.extract(r'(chr[^\s:]*)$', expand=False)

    if parsed["chromosome"].isnull().any() or parsed["start"].isnull().any() or parsed["end"].isnull().any():
        bad_examples = peak_ids[
            parsed["chromosome"].isnull() | parsed["start"].isnull() | parsed["end"].isnull()
        ].head(3).tolist()
        raise ValueError(
            "Malformed peak IDs. Expect one of: 'chr:start-end', 'chr-start-end', "
            "or BED-like 'chrom start end ...'. "
            f"Examples that failed parsing: {bad_examples}"
        )

    peak_df = pd.DataFrame({
        # "peak_id": [f"peak{i + 1}" for i in range(len(peak_ids))],
        "chromosome": parsed["chromosome"],
        "start": pd.to_numeric(parsed["start"], errors='coerce').astype(int),
        "end": pd.to_numeric(parsed["end"], errors='coerce').astype(int),
        "strand": ["."] * len(peak_ids)
    })
    
    peak_df["peak_id"] = (
        peak_df["chromosome"].astype(str) + ":" +
        peak_df["start"].astype(str) + "-" +
        peak_df["end"].astype(str)
    )
    
    return peak_df


def get_peak_length(peak_id_col: pd.Series) -> pd.Series:    
    """
    Finds the base pair lengths for a Series of genomic ranges in chr:start-end format.

    Args:
        peak_id_col (pd.Series): Series of genomic locations in chr:start-end format.

    Returns:
        pd.Series: base pair lengths of the genomic ranges proviced.
    """
    peak_col_split = peak_id_col.str.extract(r'(chr[\w]+):([0-9]+)-([0-9]+)').dropna()
    return np.abs(peak_col_split[2].astype(int) - peak_col_split[1].astype(int))


def find_genes_near_peaks(
    peak_bed: pybedtools.BedTool, 
    tss_bed: pybedtools.BedTool, 
    tss_distance_cutoff: int|float = 1e6
    ) -> pd.DataFrame:
    """
    Identify genes whose transcription start sites (TSS) are near scATAC-seq peaks.
    
    This function:
        1. Uses BedTools to find peaks that are within tss_distance_cutoff bp of each gene's TSS.
        2. Converts the BedTool result to a pandas DataFrame.
        3. Computes the absolute distance between the peak end and gene start (as a proxy for TSS distance).
        
    Args:
        peak_bed (pybedtools.BedTool):
            BedTool object with scATAC-seq peaks.
        tss_bed (pybedtools.BedTool):
            BedTool object with gene TSS locations.
        tss_distance_cutoff (int|float): 
            The maximum distance (in bp) from a TSS to consider a peak as potentially regulatory.
        
    Returns:
        peak_tss_subset_df (pandas.DataFrame): 
            A DataFrame containing columns "peak_id", "target_id", and the scaled TSS distance "TSS_dist"
            for peak–gene pairs.
    """
    
    peak_tss_overlap = peak_bed.window(tss_bed, w=tss_distance_cutoff)
    
    cols = [
        "peak_chr", "peak_start", "peak_end", "peak_id",
        "gene_chr", "gene_start", "gene_end", "gene_id"
    ]
    peak_tss_overlap_df = peak_tss_overlap.to_dataframe(
        names=cols,
        low_memory=False
    )

    # Coerce numeric cols safely & drop malformed rows
    for c in ["peak_start", "peak_end", "gene_start", "gene_end"]:
        peak_tss_overlap_df[c] = pd.to_numeric(peak_tss_overlap_df[c], errors="coerce")
    peak_tss_overlap_df = peak_tss_overlap_df.dropna(subset=["peak_start", "peak_end", "gene_start", "gene_end"]).copy()
        
    # Calculate the absolute distance in basepairs between the peak's end and gene's start.
    distances = np.abs(peak_tss_overlap_df["peak_end"].values - peak_tss_overlap_df["gene_start"].values)
    peak_tss_overlap_df["TSS_dist"] = distances
    
    # Sort by the TSS distance (lower values imply closer proximity and therefore stronger association)
    peak_tss_overlap_df = peak_tss_overlap_df.sort_values("TSS_dist")
    
    return peak_tss_overlap_df


def set_tg_as_closest_gene_tss_to_peak(tf_peak_edge_df: pd.DataFrame, peaks_gene_distance_file: str) -> pd.DataFrame:
    """
    set the target gene (TG) for each TF-peak edge as the closest gene to the peak based on TSS distance.
    
    Parameters
    ----------
    tf_peak_edge_df : pd.DataFrame
        DataFrame containing TF-peak edges with a "peak_id" column.
    peaks_gene_distance_file : str
        Path to a parquet file containing peak-to-gene distance information with columns "peak_id", "target_id", and "TSS_dist_score".
        
    Returns
    -------
    pd.DataFrame
        DataFrame with the target gene (TG) for each TF-peak edge set as the closest gene to the peak based on TSS distance.
    """
    
    assert "peak_id" in tf_peak_edge_df.columns, \
        f"'peak_id' column not in tf_peak_edge_df. Columns: {tf_peak_edge_df.columns}"
    
    # Read in the peaks to TG data and pick the closest gene for each peak (maximum TSS distance score)
    peaks_near_genes_df = pd.read_parquet(peaks_gene_distance_file, engine="pyarrow")
    
    assert "target_id" in peaks_near_genes_df.columns, \
        f"'target_id' column not in peaks_gene_distance_file DataFrame. Columns: {peaks_near_genes_df.columns}"
        
    assert "TSS_dist_score" in peaks_near_genes_df.columns, \
        f"'TSS_dist_score' column not in peaks_gene_distance_file DataFrame. Columns: {peaks_near_genes_df.columns}"

    closest_gene_to_peak_df = peaks_near_genes_df.sort_values("TSS_dist_score", ascending=False).groupby("peak_id").first()
    closest_gene_to_peak_df = closest_gene_to_peak_df[["target_id"]].reset_index()

    # Set the TG for each TF-peak edge as the closest gene to the peak
    tf_peak_tg_edge_df = pd.merge(tf_peak_edge_df, closest_gene_to_peak_df, on=["peak_id"], how="left")
    return tf_peak_tg_edge_df