import torch
import torch.nn.functional as F
import numpy as np
from Bio import AlignIO
from typing import Tuple, Optional, Dict
import gzip
from sklearn.decomposition import PCA
import json
import os
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from itertools import combinations
from tqdm import tqdm
import umap.umap_ as umap


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

def compute_weights_blockwise(
    X_idx: np.ndarray,
    theta: float | None,
    gap_idx: int | None = None,
    count_gaps_as_match: bool = False,
    block_size: int = 512,
) -> Tuple[np.ndarray, float]:
    """
    Sequence reweighting with identity threshold theta in [0,1],
    computed in memory-efficient blockwise manner.
    
    Args:
        X_idx: (M, L) integer MSA index matrix
        theta: distance threshold
        gap_idx: symbol used for gaps (or None if no gaps)
        count_gaps_as_match: whether gap-gap counts as valid & equal
        block_size: number of sequences per block to process
    Returns:
        W: per-sequence weights, shape (M,)
        M_eff: float effective number of sequences
    """
    M, L = X_idx.shape
    if theta is None or theta <= 0:
        W = np.ones(M, dtype=np.float64)
        return W, float(M)
    
    similar_counts = np.zeros(M, dtype=np.int64)
    
    for i in range(0, M, block_size):
        i_end = min(i + block_size, M)
        X_i = X_idx[i:i_end]  # shape (B, L)
        
        for j in range(0, M, block_size):
            j_end = min(j + block_size, M)
            X_j = X_idx[j:j_end]  # shape (B2, L)
            
            # Broadcasted compare: (B, B2, L)
            eq = (X_i[:, None, :] == X_j[None, :, :])
            
            if gap_idx is None:
                valid = np.ones_like(eq, dtype=bool)
            else:
                valid = (X_i[:, None, :] != gap_idx) & (X_j[None, :, :] != gap_idx)
                if count_gaps_as_match:
                    gapgap = (X_i[:, None, :] == gap_idx) & (X_j[None, :, :] == gap_idx)
                    valid |= gapgap
                    eq = np.where(gapgap, True, eq)
            
            matches = (eq & valid).sum(axis=-1)      # (B, B2)
            denom = valid.sum(axis=-1)               # (B, B2)
            
            with np.errstate(divide='ignore', invalid='ignore'):
                ident = matches / denom
                ident[denom == 0] = 0.0
                
                # Handle diagonal blocks - ensure self-comparison is 1.0
                if i == j:
                    diag_size = min(i_end - i, j_end - j)
                    np.fill_diagonal(ident[:diag_size, :diag_size], 1.0)
            
            # Accumulate counts of similar sequences
            similar_counts[i:i_end] += (ident >= theta).sum(axis=1)
    
    W = 1.0 / similar_counts.astype(np.float64)
    M_eff = float(W.sum())
    
    return W, M_eff

def compute_empirical_f1(X_idx: np.ndarray, W: np.ndarray, M_eff: float, q: int):
    """
    Compute empirical single-site frequency counts weighted by W
    Returns: (L, q) matrix
    """
    _, L = X_idx.shape
    f1 = np.zeros((L, q), dtype=float)
    for i in range(L):
        f1[i] = np.bincount(X_idx[:, i], weights=W, minlength=q)
    f1 /= M_eff
    return f1

def compute_empirical_f2(X_idx: np.ndarray, W: np.ndarray, M_eff: float, q: int):
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
    f2 /= M_eff
    return f2


def compute_empirical_f3(X_idx: np.ndarray, W: np.ndarray, M_eff: float, q: int, triplet_indices: np.ndarray):
    """
    Compute empirical three-site frequency counts weighted by W for specific triplets.
    
    Args:
        X_idx: Array of shape (M, L) containing sequence indices
        W: Array of shape (M,) containing weights
        q: Number of amino acid types
        triplet_indices: Array of shape (n_triplets, 3) with (i,j,k) indices
        
    Returns:
        f3: Array of shape (n_triplets, q, q, q) with weighted frequency counts
    """
    M, L = X_idx.shape
    n_triplets = triplet_indices.shape[0]
    f3 = np.zeros((n_triplets, q, q, q), dtype=float)
    
    for t_idx, (i, j, k) in enumerate(triplet_indices):
        ai = X_idx[:, i]
        aj = X_idx[:, j]
        ak = X_idx[:, k]
        
        # Convert triplet to single index: a*q^2 + b*q + c
        idx_triplets = ai * (q * q) + aj * q + ak
        counts = np.bincount(idx_triplets, weights=W, minlength=q*q*q).reshape(q, q, q)
        f3[t_idx] = counts
    f3 /= M_eff
    return f3

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

def save_training_history(history, save_dir, pf, version):
    """
    Save training history to a JSON file.
    """
    history_path = os.path.join(save_dir, f"training_history_{pf}_v{version}.json")
    with open(history_path, 'w') as f:
        json.dump(history, f, indent=4)
    print(f"Training history saved to {history_path}")

def load_all_histories(save_dir, pf, versions):
    """
    Load training histories for all versions of a model.
    """
    histories = {}
    for version in versions:
        history_path = os.path.join(save_dir, f"training_history_{pf}_v{version}.json")
        with open(history_path, 'r') as f:
            histories[version] = json.load(f)
    return histories

def plot_training_behavior(histories, metric="train_loss"):
    """
    Plot the training behavior for all versions.
    """
    plt.figure(figsize=(10, 6))
    for version, history in histories.items():
        plt.plot(history[metric], label=f"Version {version}")
    
    plt.xlabel("Epochs")
    plt.ylabel(metric.replace("_", " ").capitalize())
    plt.title(f"{metric.replace('_', ' ').capitalize()} Over Epochs")
    plt.legend()
    plt.grid()
    plt.savefig(f"out/{metric}")
    plt.show()


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

def umap_from_onehot(idx_mat: np.ndarray, n_components=2, n_neighbors=15, min_dist=0.1, metric="euclidean"):
    """
    Perform UMAP on one-hot encoded sequences.
    
    Args:
        idx_mat (np.ndarray): Integer matrix of shape (M, L).
        n_components (int): Number of embedding dimensions.
        n_neighbors (int): Size of local neighborhood (UMAP hyperparameter).
        min_dist (float): Minimum distance between points in low-dim space.
        metric (str): Distance metric (default: 'euclidean').
    
    Returns:
        Z (np.ndarray): UMAP embedding of shape (M, n_components).
        umap_model (umap.UMAP): Fitted UMAP object.
    """
    X = one_hot_for_pca(idx_mat)  # (M, L*q)
    reducer = umap.UMAP(
        init="random",
        n_components=n_components,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric=metric,
        random_state=42
    )
    Z = reducer.fit_transform(X)
    return Z, reducer

def compute_two_site_correlations(X_idx: np.ndarray, W: np.ndarray, M_eff: float, q: int = 21) -> np.ndarray:
    """
    Compute two-site correlations C_ij(a,b) = P_ij(a,b) - P_i(a)P_j(b).
    
    Args:
        sequences: Array of shape (M, L) containing sequence indices
        q: Number of amino acid types
        
    Returns:
        correlations: Array of shape (L, L, q, q) with connected correlations
    """
    M, L = X_idx.shape
    
    # Compute weighted single-site and two-site frequencies
    f1 = compute_empirical_f1(X_idx, W, M_eff, q)
    f2 = compute_empirical_f2(X_idx, W, M_eff, q)
    
    # Initialize correlation matrix
    correlations = np.zeros((L, L, q, q))
    
    for i in range(L):
        for j in range(i, L):
            outer = np.outer(f1[i], f1[j])  # shape (q, q)

            # Connected correlation
            correlations[i, j] = f2[i, j] - outer
            correlations[j, i] = correlations[i, j].T  # enforce symmetry
    
    return correlations

def compute_three_site_correlations(X_idx: np.ndarray, W: np.ndarray, M_eff, q: int = 21,
                                           max_triplets: Optional[int] = 1000) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute three-site connected correlations C_ijk(a,b,c) using weighted frequencies.
    
    Args:
        X_idx: Array of shape (M, L) containing sequence indices
        W: Array of shape (M,) containing weights
        q: Number of amino acid types
        max_triplets: Maximum number of triplets to compute (for efficiency)
        
    Returns:
        correlations: Array of shape (n_triplets, q, q, q) with connected correlations
        triplet_indices: Array of shape (n_triplets, 3) with (i,j,k) indices
    """
    M, L = X_idx.shape
    
    # Generate all possible triplets or sample if too many
    all_triplets = list(combinations(range(L), 3))
    if max_triplets and len(all_triplets) > max_triplets:
        np.random.seed(42)  # For reproducibility
        selected_indices = np.random.choice(len(all_triplets), max_triplets, replace=False)
        triplets = [all_triplets[i] for i in selected_indices]
    else:
        triplets = all_triplets
    
    triplet_indices = np.array(triplets)
    n_triplets = len(triplets)
    correlations = np.zeros((n_triplets, q, q, q))
    
    # Compute weighted frequencies
    f1 = compute_empirical_f1(X_idx, W, M_eff, q)
    f2 = compute_empirical_f2(X_idx, W, M_eff, q)
    f3 = compute_empirical_f3(X_idx, W, M_eff, q, triplet_indices)
    
    for t_idx, (i, j, k) in enumerate(tqdm(triplets, desc="Computing 3-site correlations")):
        for a in range(q):
            for b in range(q):
                for c in range(q):
                    # Connected correlation C_ijk(a,b,c) = P_ijk(a,b,c) - P_i(a)P_j(b)P_k(c)
                    # - C_ij(a,b)P_k(c) - C_ik(a,c)P_j(b) - C_jk(b,c)P_i(a)
                    
                    # Joint frequency P_ijk(a,b,c)
                    joint_freq = f3[t_idx, a, b, c]
                    
                    # Independent term: P_i(a)P_j(b)P_k(c)
                    independent_term = f1[i, a] * f1[j, b] * f1[k, c]
                    
                    # Two-site correlation terms
                    c_ij = f2[i, j, a, b] - f1[i, a] * f1[j, b]
                    c_ik = f2[i, k, a, c] - f1[i, a] * f1[k, c]
                    c_jk = f2[j, k, b, c] - f1[j, b] * f1[k, c]
                    
                    two_site_terms = (c_ij * f1[k, c] +
                                    c_ik * f1[j, b] +
                                    c_jk * f1[i, a])
                    
                    connected_corr = joint_freq - independent_term - two_site_terms
                    correlations[t_idx, a, b, c] = connected_corr
    
    return correlations, triplet_indices

def extract_nonzero_correlations(correlations: np.ndarray, 
                                threshold: float = 1e-10) -> np.ndarray:
    """
    Extract non-zero correlation values from correlation tensor.
    
    Args:
        correlations: Correlation tensor
        threshold: Minimum absolute value to consider non-zero
        
    Returns:
        Non-zero correlation values as 1D array
    """
    flat_corrs = correlations.flatten()
    return flat_corrs[np.abs(flat_corrs) > threshold]

def plot_correlation_comparison(data_corrs: np.ndarray, 
                              model_corrs: np.ndarray,
                              correlation_type: str = "Two-site",
                              figsize: Tuple[int, int] = (8, 6),
                              color: str = 'blue',
                              alpha: float = 0.6,
                              s: float = 1.0) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot correlation comparison between data and model samples.
    
    Args:
        data_corrs: Correlations from real data
        model_corrs: Correlations from model samples
        correlation_type: Type of correlation for labeling
        figsize: Figure size
        color: Point color
        alpha: Point transparency
        s: Point size
        
    Returns:
        Figure and axes objects
    """
    # Extract non-zero correlations
    data_flat = extract_nonzero_correlations(data_corrs)
    model_flat = extract_nonzero_correlations(model_corrs)
    
    # Ensure same length by taking minimum
    min_len = min(len(data_flat), len(model_flat))
    data_flat = data_flat[:min_len]
    model_flat = model_flat[:min_len]
    
    # Compute statistics
    correlation, p_value = pearsonr(data_flat, model_flat)
    slope = np.polyfit(data_flat, model_flat, 1)[0]
    
    # Create plot
    fig, ax = plt.subplots(figsize=figsize)
    
    # Scatter plot
    ax.scatter(data_flat, model_flat, alpha=alpha, s=s, color=color, 
              label=f'arDCA, Pearson: {correlation:.2f}\nSlope: {slope:.2f}')
    
    # Perfect correlation line
    min_val = min(np.min(data_flat), np.min(model_flat))
    max_val = max(np.max(data_flat), np.max(model_flat))
    ax.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.8, linewidth=1)
    
    # Fitted line
    fit_line = slope * data_flat + (np.mean(model_flat) - slope * np.mean(data_flat))
    sorted_idx = np.argsort(data_flat)
    ax.plot(data_flat[sorted_idx], fit_line[sorted_idx], color=color, alpha=0.8, linewidth=1)
    
    # Formatting
    ax.set_xlabel('Data', fontsize=12)
    ax.set_ylabel('Sample', fontsize=12)
    ax.set_title(f'{correlation_type} Connected Correlations', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')
    
    return fig, ax

def comprehensive_correlation_analysis(data_sequences: np.ndarray,
                                     model_sequences: np.ndarray,
                                     data_weights: np.ndarray,
                                     data_M_eff: float,
                                     model_weights: np.ndarray,
                                     model_M_eff: float,
                                     q: int = 21,
                                     max_triplets: Optional[int] = 500,
                                     figsize: Tuple[int, int] = (15, 5)) -> Dict:
    """
    Perform comprehensive correlation analysis and create comparison plots.
    
    Args:
        data_sequences: Real MSA sequences of shape (M_data, L)
        model_sequences: Model-generated sequences of shape (M_model, L)
        q: Number of amino acid types
        max_triplets: Maximum triplets for 3-site correlations
        figsize: Figure size for the combined plot
        
    Returns:
        Dictionary containing all computed correlations and statistics
    """
    print("Computing single-site frequencies...")
    data_f1 = compute_empirical_f1(data_sequences, data_weights, data_M_eff, q)
    model_f1 = compute_empirical_f1(model_sequences, model_weights, model_M_eff, q)
    
    print("Computing two-site correlations...")
    data_c2 = compute_two_site_correlations(data_sequences, data_weights, data_M_eff, q)
    model_c2 = compute_two_site_correlations(model_sequences, model_weights, model_M_eff, q)
    
    print("Computing three-site correlations...")
    data_c3, triplet_indices = compute_three_site_correlations(data_sequences, data_weights, data_M_eff, q, max_triplets)
    model_c3, _ = compute_three_site_correlations(model_sequences, model_weights, model_M_eff, q, max_triplets)
    
    
    # Create comparison plots
    fig, axes = plt.subplots(1, 3, figsize=figsize)
    
    # Single-site frequencies
    data_f1_flat = data_f1.flatten()
    model_f1_flat = model_f1.flatten()
    nonzero_mask = data_f1_flat > 1e-6  # Only non-zero frequencies
    
    corr_f1, _ = pearsonr(data_f1_flat[nonzero_mask], model_f1_flat[nonzero_mask])
    slope_f1 = np.polyfit(data_f1_flat[nonzero_mask], model_f1_flat[nonzero_mask], 1)[0]
    
    axes[0].scatter(data_f1_flat[nonzero_mask], model_f1_flat[nonzero_mask], 
                   alpha=0.6, s=2, color='blue')
    axes[0].plot([0, 1], [0, 1], 'k--', alpha=0.8, linewidth=1)
    axes[0].set_xlabel('Data')
    axes[0].set_ylabel('Sample')
    axes[0].set_title('$F_i$')
    axes[0].text(0.05, 0.95, f'arDCA, Pearson: {corr_f1:.1f}\nSlope: {slope_f1:.2f}', 
                transform=axes[0].transAxes, verticalalignment='top', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    axes[0].grid(True, alpha=0.3)
    
    # Two-site correlations
    data_c2_flat = extract_nonzero_correlations(data_c2)
    model_c2_flat = extract_nonzero_correlations(model_c2)
    min_len = min(len(data_c2_flat), len(model_c2_flat))
    
    corr_c2, _ = pearsonr(data_c2_flat[:min_len], model_c2_flat[:min_len])
    slope_c2 = np.polyfit(data_c2_flat[:min_len], model_c2_flat[:min_len], 1)[0]
    
    axes[1].scatter(data_c2_flat[:min_len], model_c2_flat[:min_len], 
                   alpha=0.4, s=1, color='blue')
    min_val = min(np.min(data_c2_flat[:min_len]), np.min(model_c2_flat[:min_len]))
    max_val = max(np.max(data_c2_flat[:min_len]), np.max(model_c2_flat[:min_len]))
    axes[1].plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.8, linewidth=1)
    axes[1].set_xlabel('Data')
    axes[1].set_ylabel('Sample')
    axes[1].set_title('$C_{ij}$')
    axes[1].text(0.05, 0.95, f'arDCA, Pearson: {corr_c2:.1f}\nSlope: {slope_c2:.2f}', 
                transform=axes[1].transAxes, verticalalignment='top', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    axes[1].grid(True, alpha=0.3)
    
    # Three-site correlations
    data_c3_flat = extract_nonzero_correlations(data_c3)
    model_c3_flat = extract_nonzero_correlations(model_c3)
    min_len = min(len(data_c3_flat), len(model_c3_flat))
    
    corr_c3, _ = pearsonr(data_c3_flat[:min_len], model_c3_flat[:min_len])
    slope_c3 = np.polyfit(data_c3_flat[:min_len], model_c3_flat[:min_len], 1)[0]
    
    axes[2].scatter(data_c3_flat[:min_len], model_c3_flat[:min_len], 
                   alpha=0.4, s=1, color='blue')
    min_val = min(np.min(data_c3_flat[:min_len]), np.min(model_c3_flat[:min_len]))
    max_val = max(np.max(data_c3_flat[:min_len]), np.max(model_c3_flat[:min_len]))
    axes[2].plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.8, linewidth=1)
    axes[2].set_xlabel('Data')
    axes[2].set_ylabel('Sample')
    axes[2].set_title('$C_{ijk}$')
    axes[2].text(0.05, 0.95, f'arDCA, Pearson: {corr_c3:.1f}\nSlope: {slope_c3:.2f}', 
                transform=axes[2].transAxes, verticalalignment='top', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    axes[2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    # Return results dictionary
    results = {
        'single_site_freqs': {
            'data': data_f1,
            'model': model_f1,
            'pearson': corr_f1,
            'slope': slope_f1
        },
        'two_site_corrs': {
            'data': data_c2,
            'model': model_c2,
            'pearson': corr_c2,
            'slope': slope_c2
        },
        'three_site_corrs': {
            'data': data_c3,
            'model': model_c3,
            'pearson': corr_c3,
            'slope': slope_c3,
            'triplet_indices': triplet_indices
        },
        'figure': fig
    }
    
    return results