import numpy as np
import pyfaidx
import tqdm
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

def onehot_dna_sequence(seq):
    """
    Fast one-hot encoding for DNA sequence.

    Returns
    -------
    np.ndarray
        Shape: (L, 4), dtype float32
    """
    
    # Fast reusable nucleotide lookup table
    _NUC_TO_IDX = np.full(256, -1, dtype=np.int16)
    _NUC_TO_IDX[ord("A")] = 0
    _NUC_TO_IDX[ord("C")] = 1
    _NUC_TO_IDX[ord("G")] = 2
    _NUC_TO_IDX[ord("T")] = 3
    _NUC_TO_IDX[ord("a")] = 0
    _NUC_TO_IDX[ord("c")] = 1
    _NUC_TO_IDX[ord("g")] = 2
    _NUC_TO_IDX[ord("t")] = 3
    
    seq_arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    idx = _NUC_TO_IDX[seq_arr]

    onehot = np.zeros((len(seq), 4), dtype=np.float32)
    valid = idx >= 0
    onehot[np.arange(len(seq))[valid], idx[valid]] = 1.0

    return onehot

def parse_peak(peak: str) -> tuple[str, int, int]:
    """
    Parse peak string like chr1:100-200.
    
    Parameters
    ----------
    peak : str
        Peak string in the format "chrom:start-end".
    
    Returns
    -------
    tuple[str, int, int]
        (chromosome, start, end)
    """
    chrom, coords = peak.split(":")
    start, end = coords.split("-")
    
    return chrom, int(start), int(end)

def load_peak_sequence(genome_fasta, selected_peak):
    """
    Load the DNA sequence for a given peak.

    Parameters
    ----------
    genome_fasta : str | Path
        Path to the genome fasta file.
    selected_peak : str
        Peak string in the format "chrom:start-end".

    Returns
    -------
    str
        The DNA sequence for the peak.
    """
    peak_chrom, peak_start, peak_end = parse_peak(selected_peak)

    # Load peak sequence using the genome fasta file
    with pyfaidx.Fasta(genome_fasta) as genome:
        peak_sequence = genome[peak_chrom][peak_start:peak_end].seq.upper()
        
    return peak_sequence

def load_chrom_sizes(chromsizes_file):
    """
    Load chromosome sizes from a chrom.sizes file.

    Parameters
    ----------
    chromsizes_file : str | Path
        Path to the chrom.sizes file.

    Returns
    -------
    dict
        Dictionary mapping chromosome names to sizes.
    """
    chrom_sizes = {}
    
    with open(chromsizes_file, "r") as f:
        for line in f:
            chrom, size_str = line.strip().split("\t")
            chrom_sizes[chrom] = int(size_str)
    
    return chrom_sizes

def _centered_peak_to_onehot(
    peak_id: str,
    genome,
    chrom_sizes: dict[str, int],
    flank_size: int,
    dtype=np.uint8,
    pad_out_of_bounds: bool = True,
):
    """
    Encode one centered peak window into a one-hot DNA matrix.

    Returns
    -------
    np.ndarray
        Shape [2 * flank_size, 4]
    """
    chrom, peak_start, peak_end = parse_peak(peak_id)

    if chrom not in chrom_sizes:
        raise KeyError(
            f"Chromosome {chrom!r} not found in chrom_sizes. "
            f"Peak: {peak_id}"
        )

    chrom_size = chrom_sizes[chrom]
    seq_len = 2 * flank_size

    peak_center = (peak_start + peak_end) // 2
    seq_start = peak_center - flank_size
    seq_end = peak_center + flank_size

    fetch_start = seq_start
    fetch_end = seq_end

    left_pad = 0
    right_pad = 0

    if fetch_start < 0:
        left_pad = -fetch_start
        fetch_start = 0

    if fetch_end > chrom_size:
        right_pad = fetch_end - chrom_size
        fetch_end = chrom_size

    if fetch_end < fetch_start:
        fetch_end = fetch_start

    seq = genome[chrom][fetch_start:fetch_end].seq.upper()

    if pad_out_of_bounds:
        if left_pad:
            seq = ("N" * left_pad) + seq
        if right_pad:
            seq = seq + ("N" * right_pad)

        if len(seq) < seq_len:
            seq = seq + ("N" * (seq_len - len(seq)))
        elif len(seq) > seq_len:
            seq = seq[:seq_len]
    else:
        if len(seq) != seq_len:
            raise ValueError(
                f"Peak {peak_id} produced sequence length {len(seq)}, "
                f"but expected {seq_len}. Use pad_out_of_bounds=True "
                f"for fixed-length output."
            )

    onehot = onehot_dna_sequence(seq).astype(dtype, copy=False)

    if onehot.shape != (seq_len, 4):
        raise ValueError(
            f"Peak {peak_id} produced one-hot shape {onehot.shape}, "
            f"but expected {(seq_len, 4)}."
        )

    return onehot


def _init_genome_handle(genome_fasta: str) -> None:
    global _GENOME_HANDLE
    _GENOME_HANDLE = pyfaidx.Fasta(genome_fasta)


def _encode_peak_chunk(args):
    """
    Worker function for multiprocessing.

    Each worker opens the FASTA once per chunk, not once per peak.
    """
    (
        peak_chunk,
        genome_fasta,
        chrom_sizes,
        flank_size,
        dtype,
        pad_out_of_bounds,
    ) = args

    results = []

    genome = _GENOME_HANDLE

    for peak_id in peak_chunk:
        onehot = _centered_peak_to_onehot(
            peak_id=peak_id,
            genome=genome,
            chrom_sizes=chrom_sizes,
            flank_size=flank_size,
            dtype=dtype,
            pad_out_of_bounds=pad_out_of_bounds,
        )
        results.append((peak_id, onehot))

    return results


def _iter_chunks(items, chunk_size: int):
    """
    Yield lists of up to chunk_size items.
    """
    chunk = []

    for item in items:
        chunk.append(item)

        if len(chunk) >= chunk_size:
            yield chunk
            chunk = []

    if chunk:
        yield chunk
        

def create_centered_peak_onehot_array(
    peak_ids: list[str],
    genome_fasta: str | Path,
    chrom_sizes: dict[str, int],
    peak_id_to_idx: dict[str, int],
    flank_size: int,
    dtype=np.uint8,
    pad_out_of_bounds: bool = True,
    show_progress: bool = True,
    num_workers: int = 1,
    chunk_size: int = 1000,
):
    """
    Create a stacked one-hot encoded DNA array using an existing peak_id_to_idx map.

    Parameters
    ----------
    peak_ids : list[str]
        Peak IDs to encode. These should all exist in peak_id_to_idx.
    genome_fasta : str | Path
        Path to genome FASTA.
    chrom_sizes : dict[str, int]
        Dictionary mapping chromosome names to chromosome sizes.
    peak_id_to_idx : dict[str, int]
        Existing mapping from peak_id -> row index.
    flank_size : int
        Number of bases on each side of the peak center.
        Output length is 2 * flank_size.
    dtype : numpy dtype
        Output dtype. np.uint8 is recommended for one-hot DNA.
    pad_out_of_bounds : bool
        Whether to pad with N if the requested window goes out of bounds.
    show_progress : bool
        Whether to show tqdm progress bar.
    num_workers : int
        Number of worker processes to use. Use 1 to run serially.
    chunk_size : int
        Number of peaks per worker task when num_workers > 1.

    Returns
    -------
    np.ndarray
        Array of shape [len(peak_id_to_idx), 2 * flank_size, 4].
    """
    genome_fasta = Path(genome_fasta)

    if not genome_fasta.exists():
        raise FileNotFoundError(f"Genome FASTA file not found: {genome_fasta}")

    if flank_size is None:
        raise ValueError("flank_size must be provided for a stacked array.")

    peak_ids = list(peak_ids)

    missing_peaks = [
        peak_id for peak_id in peak_ids
        if peak_id not in peak_id_to_idx
    ]

    if missing_peaks:
        raise KeyError(
            f"{len(missing_peaks)} peak_ids are missing from peak_id_to_idx. "
            f"Example: {missing_peaks[:5]}"
        )

    seq_len = 2 * flank_size
    num_output_peaks = len(peak_id_to_idx)
    num_encoded_peaks = len(peak_ids)

    peak_onehot_array = np.zeros(
        (num_output_peaks, seq_len, 4),
        dtype=dtype,
    )

    pbar_kwargs = dict(
        total=num_encoded_peaks,
        desc="One-hot peaks",
        disable=not show_progress,
        dynamic_ncols=True,
        miniters=max(num_encoded_peaks // 1000, 1),
    )

    if num_workers <= 1:
        with pyfaidx.Fasta(str(genome_fasta)) as genome:
            for peak_id in tqdm(peak_ids, **pbar_kwargs):
                peak_idx = peak_id_to_idx[peak_id]

                peak_onehot_array[peak_idx] = _centered_peak_to_onehot(
                    peak_id=peak_id,
                    genome=genome,
                    chrom_sizes=chrom_sizes,
                    flank_size=flank_size,
                    dtype=dtype,
                    pad_out_of_bounds=pad_out_of_bounds,
                )

    else:
        if chunk_size <= 0:
            raise ValueError("chunk_size must be > 0.")

        chunk_iter = _iter_chunks(peak_ids, chunk_size)

        task_iter = (
            (
                peak_chunk,
                genome_fasta,
                chrom_sizes,
                flank_size,
                dtype,
                pad_out_of_bounds,
            )
            for peak_chunk in chunk_iter
        )

        with ProcessPoolExecutor(
            max_workers=num_workers,
            initializer=_init_genome_handle,
            initargs=(str(genome_fasta),),
        ) as executor:
            encoded_chunk_iter = executor.map(
                _encode_peak_chunk,
                task_iter,
                chunksize=1,
            )

            with tqdm(**pbar_kwargs) as pbar:
                for encoded_chunk in encoded_chunk_iter:
                    for peak_id, onehot in encoded_chunk:
                        peak_idx = peak_id_to_idx[peak_id]
                        peak_onehot_array[peak_idx] = onehot

                    pbar.update(len(encoded_chunk))

    return peak_onehot_array
    