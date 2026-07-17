
import duckdb
import requests
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
import logging
from pathlib import Path
import numpy as np

def fetch_chip_atlas_tf(
    tf: str, 
    genome: str, 
    threshold: str = "05", 
    timeout: int = 120,
    out_dir: Path | None = None
    ) -> tuple[str, pd.DataFrame | None, Exception | None]:
    """
    Fetch ChIP-Atlas TF peaks for a given transcription factor and genome assembly.
    
    Parameters
    ----------
    tf : str
        Transcription factor name (e.g., "CTCF").
    genome : str
        Genome assembly (e.g., "mm10" for mouse, "hg38" for human).
    threshold : str
        ChIP-Atlas significance threshold (default: "05").
    timeout : int
        Timeout in seconds for the HTTP request (default: 120).
    out_dir : Path | None
        Optional directory to save the resulting DataFrame as a parquet file.
        
    Returns
    -------
    tuple
        (tf, pd.DataFrame of peaks, error if any)
    """
    
    tf_canon = tf.replace("-", "")
    
    url = (
        f"https://chip-atlas.dbcls.jp/data/{genome}/assembled/"
        f"Oth.ALL.{threshold}.{tf_canon}.AllCell.bed"
    )

    try:
        with requests.get(url, stream=True, timeout=timeout) as r:
            r.raise_for_status()

            df = pd.read_csv(
                r.raw,
                sep="\t",
                comment="t",
                header=None,
                usecols=[0, 1, 2],
                names=["peak_chr", "peak_start", "peak_end"],
                dtype={
                    "peak_chr": "category",
                    "peak_start": "int32",
                    "peak_end": "int32",
                },
            )

            if df.empty:
                return tf, None, "empty dataframe"

            # Deduplicate before writing.
            # This is much cheaper than one giant global dedup.
            df = df.drop_duplicates()

            df["source_id"] = tf

            # Keep peak coordinates separate for now.
            # Building millions of strings is expensive.
            df = df[["source_id", "peak_chr", "peak_start", "peak_end"]]

            if out_dir is not None:
                out_file = out_dir / f"{tf}.parquet"
                df.to_parquet(out_file, index=False)

            return tf, df, None

    except requests.exceptions.HTTPError as e:
        return tf, None, e
    except Exception as e:
        return tf, None, e

def fetch_chip_atlas_tf_list_to_parquet(
    tf_list: list[str],
    genome: str="mm10",
    out_dir: str="chip_atlas_tf_parquet",
    num_workers: int=10,
    threshold: str="05",
    timeout: int=120,
):
    """
    Fetch ChIP-Atlas TF peaks for a list of transcription factors and save them to parquet files.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    failed_tfs = {}
    
    existing_files = {f.stem: f for f in out_dir.glob("*.parquet")}
    if len(existing_files) > 0:
        logging.info(f"Found {len(existing_files)} / {len(tf_list)} existing parquet files. Skipping these TFs.")
    
    tf_list = [tf for tf in tf_list if tf not in existing_files]
    
    def fetch_and_save(tf):
        tf, df, error = fetch_chip_atlas_tf(
            tf=tf,
            genome=genome,
            threshold=threshold,
            timeout=timeout,
            out_dir=out_dir
        )
        return tf, df, error

    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = {
            executor.submit(fetch_and_save, tf): tf
            for tf in tf_list
        }

        for future in as_completed(futures):
            tf, df, error = future.result()

            if error is not None:
                failed_tfs[tf] = error
                logging.info(f"TF '{tf}' not found or failed: {error}")
                continue

            logging.info(f"Wrote {tf} to {out_dir / f'{tf}.parquet'}")

    return failed_tfs

def build_chip_atlas_df_from_parquet(
    parquet_dir="chip_atlas_tf_parquet",
    output_file="chip_atlas_tf_peak_edges.parquet",
):
    parquet_dir = Path(parquet_dir)

    query = f"""
    COPY (
        SELECT DISTINCT
            source_id,
            peak_chr || ':' || peak_start::VARCHAR || '-' || peak_end::VARCHAR AS peak_id
        FROM read_parquet('{parquet_dir}/*.parquet')
    )
    TO '{output_file}'
    (FORMAT PARQUET);
    """

    duckdb.sql(query)

    return output_file

def create_organism_chip_atlas_file(
    species: str, 
    ground_truth_dir: Path, 
    tf_chip_seq_save_dir: Path,
    tf_names: np.ndarray,
    num_workers: int = 10
    ) -> pd.DataFrame:
    
    full_chip_atlas_path = ground_truth_dir / f"chip_atlas_{species}_all.parquet"
    
    if not Path(full_chip_atlas_path).exists():
        fetch_chip_atlas_tf_list_to_parquet(
            tf_names, 
            genome=species, 
            out_dir=tf_chip_seq_save_dir,
            num_workers=num_workers
            )
        
        build_chip_atlas_df_from_parquet(
            parquet_dir=tf_chip_seq_save_dir, 
            output_file=full_chip_atlas_path
        )
        
        chip_atlas_full_df: pd.DataFrame = pd.read_parquet(full_chip_atlas_path)
        
        logging.info(f"Fetched {len(chip_atlas_full_df)} TF-DNA interactions from ChIP-Atlas for {species}. Saving...")
        chip_atlas_full_df.to_parquet(full_chip_atlas_path, index=False)
    else:
        chip_atlas_full_df: pd.DataFrame = pd.read_parquet(full_chip_atlas_path)
        logging.info(f"Loaded {len(chip_atlas_full_df)} TF-DNA interactions from existing file for {species}.")

    return chip_atlas_full_df