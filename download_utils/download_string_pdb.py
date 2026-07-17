    

from typing import Dict
from pathlib import Path
import os
import logging
import requests

def download_string_v12_files(string_dir: str, string_org_code: str) -> Dict[str, str]:
    """
    Download STRING v12.0 files if they don't exist locally.

    Parameters
    ----------
    string_dir : str
        Directory to store STRING files.
    string_org_code : str
        NCBI taxonomy code used by STRING
        (e.g., '10090' for mouse, '9606' for human).

    Returns
    -------
    dict
        {
            'protein_info_gz',
            'protein_links_detailed_gz',
            'protein_info_url',
            'protein_links_detailed_url',
        }
    """
    base = "https://stringdb-downloads.org/download"
    Path(string_dir).mkdir(parents=True, exist_ok=True)

    files = {
        "protein_info_gz": f"{string_org_code}.protein.info.v12.0.txt.gz",
        "protein_links_detailed_gz": f"{string_org_code}.protein.links.detailed.v12.0.txt.gz",
    }

    urls = {
        "protein_info_url": f"{base}/protein.info.v12.0/{files['protein_info_gz']}",
        "protein_links_detailed_url": f"{base}/protein.links.detailed.v12.0/{files['protein_links_detailed_gz']}",
    }

    paths = {k: os.path.join(string_dir, v) for k, v in files.items()}

    def _download(url: str, dest_path: str, chunk_size: int = 1 << 20) -> None:
        """
        Stream-download a file from `url` to `dest_path` safely.

        - Creates directories as needed
        - Writes to a temporary file first, then renames atomically
        - Logs progress
        """
        dest = Path(dest_path)
        dest.parent.mkdir(parents=True, exist_ok=True)

        tmp_path = dest.with_suffix(dest.suffix + ".tmp")
        logging.info(f"   Downloading {url} → {dest_path}")

        with requests.get(url, stream=True, timeout=60) as r:
            r.raise_for_status()
            with open(tmp_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=chunk_size):
                    if chunk:
                        f.write(chunk)

        tmp_path.replace(dest)
        logging.info(f"   Download complete: {dest.resolve()}")

    # Download if missing
    for file_key, url_key in zip(files.keys(), urls.keys()):
        path = paths[file_key]
        url = urls[url_key]
        if not os.path.exists(path):
            _download(url, path)
        else:
            logging.info(f"   Found existing: {path}")

    return {**paths, **urls}