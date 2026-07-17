from pathlib import Path
import time
import logging
from Bio import Entrez, SeqIO

def download_gene_protein_fastas(
    gene_names: list[str],
    organism: str,
    output_dir: str | Path,
    email: str,
    api_key: str | None = None,
    retmax: int = 25,
    delay: float = 0.5,
    max_tries: int = 3,
    sleep_between_tries: float = 15,
):
    """
    Download one representative RefSeq protein FASTA per gene.

    Saves:
        output_dir/{gene_name}_protein.fasta

    Uses a delay between genes to avoid NCBI rate-limit issues.
    
    Parameters
    ----------
    gene_names : list[str]
        List of gene names to download.
    organism : str
        Organism name (e.g., "Mus musculus").
    output_dir : str | Path
        Directory to save FASTA files.
    email : str
        Email address for NCBI Entrez API.
    api_key : str | None
        Optional NCBI Entrez API key for higher rate limits.
    retmax : int
        Maximum number of records to fetch per gene (default: 25).
    delay : float
        Delay in seconds between processing each gene (default: 0.5).
    max_tries : int
        Maximum number of retries for Entrez requests (default: 3).
    sleep_between_tries : float
        Seconds to sleep between retries (default: 15).
        
    Returns
    -------
    dict
        Mapping of gene_name -> Path to saved FASTA file (or None if failed).
    """

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    Entrez.email = email
    Entrez.max_tries = max_tries
    Entrez.sleep_between_tries = sleep_between_tries

    if api_key is not None:
        Entrez.api_key = api_key

    saved_files = {}
    
    # Check if the gene names are already downloaded to avoid unnecessary API calls
    available_files = {f.stem.replace("_protein", ""): f for f in output_dir.glob("*_protein.fasta")}
    gene_names = [gene for gene in gene_names if gene not in available_files]

    if not gene_names:
        logging.info("All gene FASTA files already exist. No downloads needed.")
        return {gene: available_files[gene] for gene in gene_names}
    
    for i, gene_name in enumerate(gene_names, start=1):
        search_term = (
            f'{gene_name}[Gene Name] '
            f'AND {organism}[Organism] '
            f'AND srcdb_refseq[PROP]'
        )

        try:
            with Entrez.esearch(
                db="protein",
                term=search_term,
                retmax=retmax,
            ) as search_handle:
                search_results = Entrez.read(search_handle)

            protein_ids = search_results.get("IdList", [])

            if not protein_ids:
                logging.info(f"[{i}/{len(gene_names)}] No records found for {gene_name}")
                saved_files[gene_name] = None
                time.sleep(delay)
                continue

            with Entrez.efetch(
                db="protein",
                id=protein_ids,
                rettype="gb",
                retmode="text",
            ) as fetch_handle:
                records = list(SeqIO.parse(fetch_handle, "genbank"))

            if not records:
                logging.info(f"[{i}/{len(gene_names)}] Could not parse records for {gene_name}")
                saved_files[gene_name] = None
                time.sleep(delay)
                continue

            def protein_rank(record):
                accession = record.id
                description = record.description.lower()
                keywords = [k.lower() for k in record.annotations.get("keywords", [])]

                is_refseq_select = (
                    "refseq select" in description
                    or "refseq select" in keywords
                )

                is_np = accession.startswith("NP_")
                is_xp = accession.startswith("XP_")
                is_low_quality = "low quality protein" in description

                return (
                    not is_refseq_select,
                    not is_np,
                    is_xp,
                    is_low_quality,
                    -len(record.seq),
                )

            best_record = sorted(records, key=protein_rank)[0]

            output_file = output_dir / f"{gene_name}_protein.fasta"

            with open(output_file, "w") as f:
                SeqIO.write(best_record, f, "fasta")

            saved_files[gene_name] = output_file

            logging.info(
                f"[{i}/{len(gene_names)}] Saved {gene_name}: "
                f"{best_record.id} ({len(best_record.seq)} aa)"
            )

        except Exception as e:
            logging.info(f"[{i}/{len(gene_names)}] Failed for {gene_name}: {e}")
            saved_files[gene_name] = None

        time.sleep(delay)

    return saved_files