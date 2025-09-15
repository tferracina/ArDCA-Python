import torch
import torch.nn.functional as F
import numpy as np
from Bio import AlignIO
from typing import Tuple, Union
import gzip
from sklearn.decomposition import PCA
from torch.utils.data import DataLoader, TensorDataset


# -- loading and processing data --
q = 21

ALPHABET = "-ACDEFGHIKLMNPQRSTVWY"
AA2IDX = {aa:i for i, aa in enumerate(ALPHABET)}

def aa2idx(aa: str) -> int:
    return AA2IDX.get(aa.upper(), 0)

def read_fasta_alignment(filename: str, max_gap_fraction: float,  max_col_gap_fraction=None) -> np.ndarray:
    """Parses an MSA in FASTA file -> returns matrix of intergers 
    discard if seq contains a fraction of gaps >  `max_gap_fraction` - set 1 to keep all
    discard if column contains a fraction of gaps > `max_col_gap_fraction` - if set
    """
    if filename.endswith('.gz'):
        with gzip.open(filename, 'rt') as handle:
            alignment = AlignIO.read(handle, "fasta")
    else:
        alignment = AlignIO.read(filename, "fasta")

    L = alignment.get_alignment_length()

    # Filter sequences based on gap fraction
    keep = []
    for rec in alignment:
        s = str(rec.seq).upper()
        if s.count('-') / L <= max_gap_fraction:
            keep.append(s)
    if not keep:
        raise ValueError("No sequences passed gap filter (max_gap_fraction)={max_gap_fraction})")
    
    # Create numerical matrix
    M = len(keep)
    idx_matrix = np.zeros((M, L), dtype = np.int16)

    for i, s in enumerate(keep):
        idx_matrix[i] = [aa2idx(c) for c in s]

    if max_col_gap_fraction is not None:
        col_gap_frac = (idx_matrix == 0).mean(axis=0)
        col_mask = col_gap_frac <= max_col_gap_fraction
        idx_matrix = idx_matrix[:, col_mask]

    return idx_matrix

def encode_sequence(X_idx: np.ndarray, q: int = 21, device: str = "cpu") -> torch.Tensor:
    """
    Convert index-encoded sequences [M,L] into one-hot [M,L,q].
    """
    A = torch.as_tensor(X_idx, dtype=torch.long, device=device)
    X_one_hot = F.one_hot(A, num_classes=q).to(torch.float32)
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

    # gap-gap comparisons
    if gap_idx is None:
        valid = np.ones_like(eq, dtype=bool)
    else:
        valid = (X_idx[:, None, :] != gap_idx) & (X_idx[None, :, :] != gap_idx)
        if count_gaps_as_match:
            # treat gap-gap as valid & equal
            gapgap = (X_idx[:, None, :] == gap_idx) & (X_idx[None, :, :] == gap_idx)
            valid = valid | gapgap
            eq = np.where(gapgap, True, eq)

    matches = (eq & valid).sum(axis=-1) # equality
    denom = valid.sum(axis=-1) # number of valid positions 
    with np.errstate(divide='ignore', invalid='ignore'):
        ident = matches / denom # check seqid
        ident[denom == 0] = 0.0  # fallback for overlapping nongap sites

    np.fill_diagonal(ident, 1.0)                
    similar_counts = (ident >= theta).sum(axis=1) 
    W = 1.0 / similar_counts
    M_eff = float(W.sum())
    return W, M_eff

def compute_empirical_f1(X_idx: np.ndarray, W: np.ndarray, q: int):
    """
    Compute empirical single-site frequency counts weighted by W
    Returns: (L, q) matrix
    """
    _, L = X_idx.shape
    f1 = np.zeros((L, q), dtype=float)
    for i in range(L):
        f1[i] = np.bincount(X_idx[:, i], weights=W, minlength=q)
    return f1

def compute_empirical_f2(X_idx: np.ndarray, W: np.ndarray, q: int):
    """
    Compute empirical pairwise frequency counts weighted by W
    Returns: (L, L, q, q) matrix
    """
    _, L = X_idx.shape
    f2 = np.zeros((L, L, q, q), dtype=float)

    for i in range(L):
        ai = X_idx[:, i]
        for j in range(i, L):
            aj = X_idx[:, j]
            idx_pairs = ai * q + aj
            counts = np.bincount(idx_pairs, weights=W, minlength=q*q).reshape(q, q)
            f2[i, j] = counts
            if j != i:
                f2[j, i] = counts.T

    return f2


# -- model training utilities --
def split_sequences(X: np.ndarray, W: np.ndarray, val_frac: float = 0.2, seed: int = 0):
    """
    Random split of sequences + weights into train/val.
    Returns (X_train, W_train, X_val, W_val)
    """
    rng = np.random.default_rng(seed)
    M = X.shape[0]
    idx = rng.permutation(M)
    m_val = max(1, int(round(val_frac * M)))
    val_idx, train_idx = idx[:m_val], idx[m_val:]
    return X[train_idx], W[train_idx], X[val_idx], W[val_idx]


# -- analyzing results --
def one_hot_for_pca(idx_mat: np.ndarray, q: int = 21) -> np.ndarray:
    M, L = idx_mat.shape
    onehot = np.eye(q, dtype=np.float32)[idx_mat]      # (M, L, q)
    X = onehot.reshape(M, L * q)                       # (M, L*q)
    return X

def pca_from_onehot(idx_mat: np.ndarray, n_components=2):
    X = one_hot_for_pca(idx_mat)            # (M, L*q)
    X -= X.mean(axis=0, keepdims=True)
    pca = PCA(n_components=n_components, svd_solver="auto")
    Z = pca.fit_transform(X)                   # (M, n_components)
    return Z, pca