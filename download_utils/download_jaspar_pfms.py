import os
from pathlib import Path
import logging
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed

def download_jaspar_pfms(save_dir: str, tax_id: str = "10090", version: int = 2024, max_workers: int = 8):
    """
    Download all JASPAR PFMs for a given organism (e.g., mouse) via REST API.

    Parameters
    ----------
    save_dir : str
        Directory to save .pfm files.
    tax_id : str
        NCBI taxonomy ID (e.g., '10090' for mouse, '9606' for human).
    version : int
        JASPAR release version.
    max_workers : int
        Parallel download threads.
    """
    os.makedirs(save_dir, exist_ok=True)

    # List endpoint for this organism
    list_url = f"https://jaspar.elixir.no/api/v1/matrix/?tax_id={tax_id}&version={version}&page_size=500"
    logging.info(f"Fetching JASPAR matrices: {list_url}")

    resp = requests.get(list_url)
    resp.raise_for_status()
    data = resp.json()

    results = data.get("results", [])
    logging.info(f"Found {len(results)} motifs for tax_id={tax_id}")

    # Build URLs for PFMs
    pfm_urls = {
        r["matrix_id"]: f"{r['url']}?format=pfm"
        for r in results if "url" in r
    }

    logging.info(f"Preparing to download {len(pfm_urls)} PFMs...")

    def _download(url: str, dest: str, chunk_size: int = 1 << 18):
        """Stream a file to disk safely."""
        dest = Path(dest)
        dest.parent.mkdir(parents=True, exist_ok=True)
        tmp = dest.with_suffix(dest.suffix + ".tmp")
        with requests.get(url, stream=True, timeout=60) as r:
            r.raise_for_status()
            with open(tmp, "wb") as f:
                for chunk in r.iter_content(chunk_size=chunk_size):
                    if chunk:
                        f.write(chunk)
        tmp.replace(dest)
        logging.info(f"Downloaded {dest.name}")

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for matrix_id, pfm_url in pfm_urls.items():
            dest = os.path.join(save_dir, f"{matrix_id}.pfm")
            if not os.path.exists(dest):
                futures.append(executor.submit(_download, pfm_url, dest))
            else:
                logging.info(f"Already exists: {dest}")
        for fut in as_completed(futures):
            fut.result()

    logging.info(f"All PFMs saved under {save_dir}")