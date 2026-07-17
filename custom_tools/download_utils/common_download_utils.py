from requests import get as request_get
from pathlib import Path
from tqdm import tqdm

def stream_download(url: str, dest: Path, chunk: int = 1 << 20, desc: str|None = None):
    """
    stream-download a file from `url` to `dest` safely, with progress bar.
    
    Parameters
    ----------
    url : str
        The URL to download from.
    dest : Path
        The destination file path.
    chunk : int
        Chunk size in bytes for streaming (default: 1MB).
    desc : Optional[str]
        Description for the progress bar. If None, uses dest.name.
    """
    
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    with request_get(url, stream=True, timeout=(10, 600)) as r:
        r.raise_for_status()
        total = int(r.headers.get("Content-Length", 0)) or None
        with open(tmp, "wb") as f, tqdm(
            total=total, unit="B", unit_scale=True, unit_divisor=1024,
            desc=desc or dest.name, dynamic_ncols=True, mininterval=0.1
        ) as pbar:
            for chunk_bytes in r.iter_content(chunk_size=chunk):
                if not chunk_bytes:
                    continue
                f.write(chunk_bytes)
                pbar.update(len(chunk_bytes))
    tmp.replace(dest)
