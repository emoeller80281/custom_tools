import pandas as pd
import subprocess
from pathlib import Path
import logging
import shutil
from typing import Union, Optional
import requests
import gzip
import pysam
from tqdm.auto import tqdm
from pybiomart import Dataset

from .common_download_utils import stream_download

def download_gene_tss_file(
    save_file: Path|str, 
    gene_dataset_name: str = "hsapiens_gene_ensembl", 
    ensembl_version: str = "useast.ensembl.org"
    ) -> pd.DataFrame:
    """
    Downloads the gene TSS coordinates from Ensembl BioMart using pybiomart.

    Parameters
    ----------
    save_file: Path|str
        Path to save the gene TSS bed file.
        
    ensembl_version : str, optional
        Ensembl host mirror to query (default: "useast.ensembl.org").
        Examples: "www.ensembl.org", "uswest.ensembl.org", etc.

    Returns
    -------
    pd.DataFrame
        A DataFrame with columns: ['gene_name', 'chromosome', 'tss']
    """

    # Connect to the Ensembl BioMart gene dataset for the organism
    dataset = Dataset(name=gene_dataset_name, host=f"http://{ensembl_version}")

    # Retrieve TSS, gene name, and chromosome
    df = dataset.query(
        attributes=[
            "external_gene_name",
            "chromosome_name",
            "transcription_start_site",
        ]
    )

    # Clean up and rename
    df = df.rename(
        columns={
            "Gene name": "name",
            "Chromosome/scaffold name": "chrom",
            "Transcription start site (TSS)": "start",
        }
    )

    # Filter out non-standard chromosomes
    df = df[df["chrom"].str.match(r"^\d+$|^X$|^Y$")]

    # Add the chr prefix to the chromosomes
    df["chrom"] = df["chrom"].astype(str)
    df["chrom"] = df["chrom"].apply(lambda c: f"chr{c}" if not c.startswith("chr") else c)

    # Drop duplicates
    df = df.drop_duplicates(subset=["name", "chrom", "start"]).reset_index(drop=True)

    df["end"] = df["start"] + 1

    df = df.dropna()


    df = df[["chrom", "start", "end", "name"]]

    save_file.parent.mkdir(parents=True, exist_ok=True)
    
    df.to_csv(save_file, sep="\t", header=False, index=False)

    return df

def download_genome_fasta(organism_code: str, save_dir: Union[str, Path]) -> Path:
    """
    Download a UCSC genome FASTA, overwrite gzip with BGZF (still .gz), and index via pysam.faidx.
    """

    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    gz_path  = save_dir / f"{organism_code}.fa.gz"      # final file (BGZF) keeps .gz suffix
    fai_path = save_dir / f"{organism_code}.fa.gz.fai"  # index alongside .gz
    gzi_path = save_dir / f"{organism_code}.fa.gz.gzi"  # BGZF index
    
    def _download_with_progress(url: str, dest: Path, chunk_size: int = 256 * 1024, desc: str = "") -> None:
        dest.parent.mkdir(parents=True, exist_ok=True)
        tmp = dest.with_suffix(dest.suffix + ".tmp")
        with requests.get(url, stream=True, timeout=(10, 600)) as r:
            r.raise_for_status()
            total = int(r.headers.get("Content-Length", 0)) or None
            with open(tmp, "wb") as f, tqdm(
                total=total,
                unit="B", unit_scale=True, unit_divisor=1024,
                desc=desc or dest.name,
                dynamic_ncols=True, mininterval=0.1,
            ) as pbar:
                for chunk in r.iter_content(chunk_size=chunk_size):
                    if not chunk:
                        continue
                    f.write(chunk)
                    pbar.update(len(chunk))
        tmp.replace(dest)

    # 1) Download gzip if missing
    if not gz_path.exists():
        url = f"https://hgdownload.soe.ucsc.edu/goldenPath/{organism_code}/bigZips/{organism_code}.fa.gz"
        logging.info(f"Downloading {organism_code} genome from:\n  {url}")
        _download_with_progress(url, gz_path, desc=gz_path.name)
        logging.info(f"  - Download complete: {gz_path}")

        
    def _is_bgzf(path: Path) -> bool:
        """
        True iff file is BGZF (gzip with 'BC' extra subfield).
        Fast path: if a .gzi sibling exists (non-empty), assume BGZF.
        """
        try:
            gzi = Path(str(path) + ".gzi")  # e.g., mm10.fa.gz.gzi
            if gzi.exists() and gzi.stat().st_size > 0:
                return True

            with open(path, "rb") as fh:
                hdr = fh.read(10)  # ID1 ID2 CM FLG MTIME(4) XFL OS  -> 10 bytes
                if len(hdr) < 10 or hdr[0:2] != b"\x1f\x8b" or hdr[2] != 8:
                    return False
                flg = hdr[3]
                if not (flg & 0x04):  # FEXTRA not set
                    return False

                # Read XLEN (2 bytes, little-endian), then that many bytes of extra
                xlen_b = fh.read(2)
                if len(xlen_b) != 2:
                    return False
                xlen = int.from_bytes(xlen_b, "little")
                extra = fh.read(xlen)

                # Walk subfields looking for 'BC'
                i = 0
                while i + 4 <= len(extra):
                    si = extra[i:i+2]
                    slen = int.from_bytes(extra[i+2:i+4], "little")
                    i += 4
                    if si == b"BC":
                        return True
                    i += slen
                return False
        except Exception:
            return False

    # 2) Ensure it's BGZF; UCSC provides plain gzip, so this will usually run once
    if not _is_bgzf(gz_path):
        logging.info(f"{gz_path.name} is plain gzip; converting to BGZF (keeping .gz)…")
        tmp_bgzf = gz_path.with_suffix(gz_path.suffix + ".bgzf.tmp")
        tmp_fa = gz_path.with_suffix("")  # uncompressed intermediate
        
        # Try using bgzip CLI (much faster, parallel) if available
        try:
            # Decompress to temp file
            logging.info(f"  - Decompressing with gunzip...")
            subprocess.run(["gunzip", "-c", str(gz_path)], 
                          stdout=open(tmp_fa, "wb"), 
                          check=True, 
                          stderr=subprocess.DEVNULL)
            
            # Compress with bgzip (parallel, much faster)
            logging.info(f"  - Compressing with bgzip...")
            subprocess.run(["bgzip", "-f", "-@", "8", str(tmp_fa)], 
                          check=True, 
                          stderr=subprocess.DEVNULL)
            
            # bgzip creates .gz file, rename to our temp name then to final
            (tmp_fa.parent / f"{tmp_fa.name}.gz").replace(tmp_bgzf)
            tmp_bgzf.replace(gz_path)
            logging.info(f"  - Converted to BGZF using bgzip")
            
        except (subprocess.CalledProcessError, FileNotFoundError):
            # Fall back to Python if bgzip not available
            logging.info(f"  - bgzip not available, using slower Python conversion...")
            if tmp_fa.exists():
                tmp_fa.unlink()
            
            # Use larger chunks (4MB) for better performance
            with gzip.open(gz_path, "rb") as fin, pysam.BGZFile(str(tmp_bgzf), "wb") as bout:
                chunk_size = 1 << 22  # 4MB chunks
                file_size = gz_path.stat().st_size
                bytes_read = 0
                last_log_pct = -1
                
                for chunk in iter(lambda: fin.read(chunk_size), b""):
                    if not chunk:
                        break
                    bout.write(chunk)
                    bytes_read += len(chunk)
                    
                    # Log progress every 10%
                    pct = int((bytes_read / file_size) * 100)
                    if pct >= last_log_pct + 10:
                        logging.info(f"    {pct}% complete")
                        last_log_pct = pct
            
            tmp_bgzf.replace(gz_path)
            logging.info(f"  - Conversion complete")

    # 3) Index the BGZF FASTA (.fai + .gzi)
    if fai_path.exists() and gzi_path.exists():
        return gz_path
    else:
        logging.info(f"Indexing {gz_path.name} with pysam.faidx …")
        pysam.faidx(str(gz_path))
        logging.info(f"  - Index created: {fai_path.name} (and {gzi_path.name})")

    return gz_path

def download_chrom_sizes(organism_code: str, save_dir: Union[str, Path]) -> Path:
    """
    Download UCSC chrom.sizes for an assembly (mm10 or hg38).

    Returns
    -------
    Path
        Path to the chrom.sizes file.
    """
    assert organism_code in ("mm10", "hg38"), \
        f"Organism code '{organism_code}' not supported (valid: 'mm10', 'hg38')."

    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    url = f"https://hgdownload.soe.ucsc.edu/goldenPath/{organism_code}/bigZips/{organism_code}.chrom.sizes"
    out_path = save_dir / f"{organism_code}.chrom.sizes"

    if out_path.exists():
        return out_path

    logging.info(f"Downloading chrom.sizes:\n  {url}")
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        tmp = out_path.with_suffix(out_path.suffix + ".tmp")
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 20):
                if chunk:
                    f.write(chunk)
        tmp.replace(out_path)

    logging.info(f"chrom.sizes saved: {out_path}")
    return out_path

def download_ncbi_gene_info(organism_code: str, out_path: Optional[Union[str, Path]] = None) -> Path:
    """
    Download NCBI Gene Info for mouse or human.
    
    Parameters
    ----------
    organism_code : str
        'mm10' for mouse or 'hg38' for human
    out_path : Optional[Union[str, Path]]
        Custom output path. If None, uses default based on organism.
    
    Returns
    -------
    Path
        Path to the downloaded gene_info.gz file
    """
    # Map organism code to NCBI species and paths
    org_map = {
        "mm10": {
            "species": "Mus_musculus",
            "subdir": "Mammalia",
            "url": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz"
        },
        "hg38": {
            "species": "Homo_sapiens",
            "subdir": "Mammalia",
            "url": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
        }
    }
    
    if organism_code not in org_map:
        raise ValueError(f"Unsupported organism_code: {organism_code}. Use 'mm10' or 'hg38'.")
    
    info = org_map[organism_code]
    url = info["url"]
    filename = f"{info['species']}.gene_info.gz"
    
    # Construct organism-specific path
    if out_path is not None:
        dest = Path(out_path)
    else:
        # Construct path relative to current working directory
        base = Path("data/genome_data/genome_annotation") / organism_code
        dest = base / filename
        # Make it absolute
        if not dest.is_absolute():
            dest = Path.cwd() / dest

    if dest.exists():
        return dest

    logging.info(f"Downloading NCBI Gene Info for {info['species']}:\n  {url}")
    stream_download(url, dest, desc=dest.name)
    logging.info(f"Saved: {dest}")
    return dest

def download_ensembl_gtf(
    organism_code: str,
    release: Optional[int] = None,
    assembly: Optional[str] = None,
    out_dir: Optional[Union[str, Path]] = None,
    decompress: bool = False,
) -> Path:
    """
    Download Ensembl GTF for mouse or human. By default keeps .gtf.gz.
    If decompress=True, also writes an uncompressed .gtf alongside.
    
    Parameters
    ----------
    organism_code : str
        'mm10' for mouse or 'hg38' for human
    release : Optional[int]
        Ensembl release version. If None, uses defaults (115 for mouse, 113 for human)
    assembly : Optional[str]
        Genome assembly. If None, uses defaults (GRCm39 for mouse, GRCh38 for human)
    out_dir : Optional[Union[str, Path]]
        Output directory. If None, uses GTF_FILE_DIR or creates default path
    decompress : bool
        Whether to also create an uncompressed .gtf file
    
    Returns
    -------
    Path
        Path to the GTF file (.gtf.gz or .gtf if decompress=True)
    """
    # Map organism code to Ensembl parameters
    org_map = {
        "mm10": {
            "species_dir": "mus_musculus",
            "species_name": "Mus_musculus",
            "default_assembly": "GRCm39",
            "default_release": 115
        },
        "hg38": {
            "species_dir": "homo_sapiens",
            "species_name": "Homo_sapiens",
            "default_assembly": "GRCh38",
            "default_release": 113
        }
    }
    
    if organism_code not in org_map:
        raise ValueError(f"Unsupported organism_code: {organism_code}. Use 'mm10' or 'hg38'.")
    
    info = org_map[organism_code]
    release = release if release is not None else info["default_release"]
    assembly = assembly if assembly is not None else info["default_assembly"]
    
    org = info["species_dir"]
    species_name = info["species_name"]
    fn_gz = f"{species_name}.{assembly}.{release}.gtf.gz"
    url = f"https://ftp.ensembl.org/pub/release-{release}/gtf/{org}/{fn_gz}"

    # Construct organism-specific path
    if out_dir is None:
        # Construct path relative to current working directory
        out_dir = Path("data/genome_data/genome_annotation") / organism_code
        # Make it absolute
        if not out_dir.is_absolute():
            out_dir = Path.cwd() / out_dir
    else:
        out_dir = Path(out_dir)
    
    out_dir.mkdir(parents=True, exist_ok=True)
    dest_gz = out_dir / fn_gz
    dest_gtf = dest_gz.with_suffix("")  # drop .gz -> .gtf

    if not dest_gz.exists():
        logging.info(f"Downloading Ensembl GTF for {species_name}:\n  {url}")
        stream_download(url, dest_gz, desc=dest_gz.name)
        logging.info(f"Saved: {dest_gz}")

    if decompress:
        if dest_gtf.exists():
            logging.info(f"Uncompressed GTF already exists: {dest_gtf}")
        else:
            logging.info(f"Decompressing → {dest_gtf.name}")
            with gzip.open(dest_gz, "rb") as fin, open(dest_gtf, "wb") as fout:
                shutil.copyfileobj(fin, fout, length=1 << 20)
            logging.info(f"Wrote: {dest_gtf}")

    return dest_gtf if decompress else dest_gz