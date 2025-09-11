import torch
import torch.nn.functional as F
import numpy as np
from Bio import AlignIO
from typing import Tuple, Union
import gzip


q = 21


def aa2idx(aa: str) -> int:
    """Return numer corresponding to amino acid, return 0 for gaps"""
    aa2num_dict = {
        'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5,
        'G': 6, 'H': 7, 'I': 8, 'K': 9, 'L': 10,
        'M': 11, 'N': 12, 'P': 13, 'Q': 14, 'R': 15,
        'S': 16, 'T': 17, 'V': 18, 'W': 19, 'Y': 20,
    }
    return aa2num_dict.get(aa, 0)


def read_fasta_alignment(filename: str, max_gap_fraction: float) -> np.ndarray:
    """Parses a FASTA file containing a MSA, and returns a matrix of intergers 

    If a seq contains a fraction of gaps that exceeds `max_gap_fraction`, it is discarded. 
    Set this value to 1 to keep all sequences.
    
    Handles both plain FASTA files and gzipped FASTA files (.gz extension).
    """
    
    # Check if file is gzipped and open accordingly
    if filename.endswith('.gz'):
        with gzip.open(filename, 'rt') as handle:
            alignment = AlignIO.read(handle, "fasta")
    else:
        alignment = AlignIO.read(filename, "fasta")

    # Filter sequences based on gap fraction
    filtered_seqs = []
    for record in alignment:
        seq_str = str(record.seq).upper()
        gap_count = seq_str.count('-')
        if gap_count / len(seq_str) <= max_gap_fraction:
            filtered_seqs.append(record)

    if not filtered_seqs:
        raise ValueError("No sequences passed gap filter (max_gap_fraction)={max_gap_fraction})")
    # Create numerical matrix
    M = len(filtered_seqs)
    L = alignment.get_alignment_length()
    matrix = np.zeros((M, L), dtype = np.int8)

    # Fill matrix
    for seq_i, record in enumerate(filtered_seqs):
        seq_str = str(record.seq).upper()
        for pos_j, char in enumerate(seq_str):
            matrix[seq_i, pos_j] = aa2idx(char)

    return matrix


def encode_sequence(X_idx: np.ndarray) -> torch.Tensor:
    """Builds X_onehot from X_idx of shape [M, L, q] (0/1 per residue)"""
    A = torch.tensor(X_idx)
    X_one_hot = F.one_hot(A.long(), num_classes=21)

    return X_one_hot


def compute_weights(X_idx: np.ndarray, 
                    theta: float | None,
                    gap_idx: int | None = None,
                    count_gaps_as_match: bool = False) -> Tuple[np.ndarray, float]:
    """
    Sequence reweighting with identity threshold theta in [0,1]
    w_m = 1 / # (seqid(m,m') >= theta) - no reweighting if theta = None
    X_idx: MSA index matrix [M, L]
    theta: Distance threshold
    """
    M, L = X_idx.shape
    if theta is None or theta <= 0:
        W = np.ones(M, dtype=np.float64)
        return W, float(M)
    
    eq = (X_idx[:, None, :] == X_idx[None, :, :]) # [M, M, L]

    # remove (or keep) gap-gap comparisons
    if gap_idx is None:
        valid = np.ones_like(eq, dtype=bool)
    else:
        valid = (X_idx[:, None, :] != gap_idx) & (X_idx[None, :, :] != gap_idx)
        if count_gaps_as_match:
            # treat gap-gap as valid & equal
            gapgap = (X_idx[:, None, :] == gap_idx) & (X_idx[None, :, :] == gap_idx)
            valid = valid | gapgap
            eq = np.where(gapgap, True, eq)

    matches = (eq & valid).sum(axis=-1) # check equality
    denom = valid.sum(axis=-1) # check number of valid positions 
    with np.errstate(divide='ignore', invalid='ignore'):
        ident = matches / denom # check sequence identity
        ident[denom == 0] = 0.0  # fallback for overlapping nongap sites

    np.fill_diagonal(ident, 1.0)                  # include self
    similar_counts = (ident >= theta).sum(axis=1) # n_m
    W = 1.0 / similar_counts
    M_eff = float(W.sum())
    return W, M_eff


def compute_empirical_f1(X_idx: np.ndarray, W: np.ndarray, q: int):
    """
    Compute empirical single-site frequency counts weighted by W
    Returns: (L, q) matrix
    """
    _, L = X_idx.shape

    f1 = np.zeros((L, q))

    for i in range(L):
        np.add.at(f1[i, :], X_idx[:, i], W)

    return f1


def compute_empirical_f2(X_idx: np.ndarray, W: np.ndarray, q: int):
    """
    Compute empirical pairwise frequency counts weighted by W
    Returns: (L, L, q, q) matrix
    """
    _, L = X_idx.shape
    f2 = np.zeros((L, L, q, q))

    for i in range(L):
        for j in range(i, L):
            idx_pairs = X_idx[:, i] * q + X_idx[:, j]
            counts = np.bincount(idx_pairs, weights=W, minlength=q*q).reshape(q, q)
            f2[i, j] = counts
            if j != i:
                f2[j, i] = counts.T

    return f2

