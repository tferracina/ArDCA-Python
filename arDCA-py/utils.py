import numpy as np
from Bio import AlignIO
from typing import Any, Tuple, Union, List
import numpy.random as random

from ar_types import ArVar, ArNet

ALL_PERMORDER = ["NATURAL", "ENTROPIC", "REV_ENTROPIC", "RANDOM"]

def aa_to_num(aa: str) -> np.int8:
    """String to number, return 21 for gaps and unrecognized capital letters"""
    aa_to_num_dict = {
        'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5,
        'G': 6, 'H': 7, 'I': 8, 'K': 9, 'L': 10,
        'M': 11, 'N': 12, 'P': 13, 'Q': 14, 'R': 15,
        'S': 16, 'T': 17, 'V': 18, 'W': 19, 'Y': 20,
    }
    return aa_to_num_dict.get(aa, 21)


def read_fasta_alignment(filename: str, max_gap_fraction: float) -> np.ndarray:
    """Parses a FASTA file containing a MSA, and returns a matrix of intergers 
    that represents one sequence per column.

    If a seq contains a fraction of gaps that exceeds `max_gap_fraction`, it is discarded. 
    Set this value to 1 to keep all sequences.
    """

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
    # Create numerical matrix (sequences as columns)
    Ngth = alignment.get_alignment_length()
    num_seqs = len(filtered_seqs)
    matrix = np.zeros((Ngth, num_seqs), dtype = np.int8)

    # Fill matrix
    for col_i, record in enumerate(filtered_seqs):
        seq_str = str(record.seq).upper()
        for row_j, char in enumerate(seq_str):
            matrix[row_j, col_i] = aa_to_num(char)

    return matrix


def remove_duplicate_sequences(matrix: np.ndarray) -> np.ndarray:
    """Remove duplicate columns"""
    matrixt = matrix.T

    _, unique_indices = np.unique(matrixt.view('i1').reshape(matrixt.shape[0], -1),
                                  axis=0, return_index=True)
    unique_indices.sort()

    return matrix[:, unique_indices]

def hamming_distance_matrix(Z: np.ndarray) -> np.ndarray:
    return np.sum(Z[:, None, :] != Z[None, :, :], axis=-1)


def compute_weights(Z: np.ndarray, max_val: int = None, theta: Union[float, str] = None) -> Tuple[np.ndarray, float]:
    """Compute the reweighting vector. Retruns vector and its sum, # of "effective" sequence
    Z: MSA Matrix (NxM)
    max_val: Maximum value in the alphabet
    theta: Distance threshold
    """
    N, M = Z.shape
    

    if q is None:
        q = Z.max()

    if not isinstance(theta, (int, float)):
        raise ValueError("theta must be int/float")

    if theta == 0:
        return np.ones(M), float(M)

    # Compute distance and weights
    thresh = np.floor(theta * N)
    distances = hamming_distance_matrix(Z.T)
    similar_counts = np.sum(distances < thresh, axis=1)
    W = 1.0 / similar_counts
    Meff = np.sum(W)

    return W, Meff


def read_fasta(filename: str, max_gap_fraction: float, theta: Any, remove_dups: bool): #add theta: Any
    """Read FASTA file"""
    Z = read_fasta_alignment(filename, max_gap_fraction)
    if remove_dups:
        Z = remove_duplicate_sequences(Z)
    N, M = Z.shape # (N length of sequence, M # of sequences)
    q = int(np.max(Z))
    W, Meff = compute_weights(Z, theta=theta)
    return W, Z, N, M, q


def checkpermorder(permorder):
    if isinstance(permorder, str):
        if not permorder in ALL_PERMORDER:
            raise ValueError("permorder not implemented")
    elif isinstance(permorder, list):
        if not is_permutation(permorder):
            raise ValueError("permorder is not a permutation")
    else:
        raise TypeError("permorder must be a string or a list of integers")


def is_permutation(permorder):
    n = len(permorder)
    seen = set()
    for x in permorder:
        if not isinstance(x, int) or x < 1 or x > n or x in seen:
            return False
        seen.add(x)
    return True


#laplace smoothing
def computep0(var: ArVar) -> np.ndarray: # MAKE SURE THAT Z IS ALSO 0-indexed !!!!!
    """
    Compute a probability distribution by accumulating weights according to category assignments,
    then applying pseudocount smoothing.
      
    Returns:
      A 1D numpy array of length q representing the smoothed probability distribution.
    """
    W, Z, q, pc = var.W, var.Z, var.q, var.pc
    p0 = np.zeros(q)
    np.add.at(p0, Z[0, :], W)

    return p0 * (1 - pc) + pc / q


def compute_empirical_freqs(Z: np.ndarray, W: np.ndarray, q: int): # MAKE SURE THAT Z IS ALSO 0-INDEXED !!!!!
    """
    Compute empirical frequency counts weighted by W for a categorical dataset.
    
     Returns:
      A (q, Nrow) matrix."""
    Nrow, _ = Z.shape

    f = np.zeros((q, Nrow))

    for i in range(Nrow):
        np.add.at(f[:, i], Z[i, :], W)

    return f


def entropy(Z: np.ndarray, W: np.ndarray):
    q = int(Z.max() + 1)
    f = compute_empirical_freqs(Z, W, q)

    epsilon = 1e-10
    mask = f > 0
    
    return -np.sum(np.where(mask, f * np.log(f + epsilon), 0), axis=0)


def unpack_params(theta: np.ndarray, var: Any) -> Tuple[np.ndarray, List[np.ndarray], List[np.ndarray]]:
    N, q = var.N, var.q

    arrJ = []  #each array corresponds to a site, holds interactiosn between pairs
    arrH = []  #each array corresponds to bias or weight for each alphabet element at a site
    counter = 0
    for site in range(N):
        _arrJ = np.zeros((q, q, site))
        for i in range(site):
            for a in range(q):
                for b in range(q):
                    counter += 1
                    _arrJ[b, a, i] = theta[counter]
        arrJ.append(_arrJ)

        _arrH = np.zeros(q)
        for a in range(q):
            counter += 1
            _arrH[a] = theta[counter]
        arrH.append(_arrH)

    assert counter == len(theta), f"Expected {len(theta)} values, but found {counter}"
    return computep0(var), arrJ, arrH


def softmax(x: np.ndarray) -> np.ndarray: #vectorized operations
    u = np.max(x, axis=0, keepdims=True)
    r = np.exp(x - u)
    r /= np.sum(r, axis=0, keepdims=True)
    return r

def softmax_inplace(x: np.ndarray):
    r = np.exp(x - np.max(x))
    r /= np.sum(r)
    return r


def sample(arnet: ArNet, msamples: int) -> np.ndarray:
    H, J, p0, idxperm = arnet.H, arnet.J, arnet.p0, arnet.idxperm
    q = p0.size
    N, _ = H.shape # where N is actually N-1 
    res_N = N + 1

    back_order = np.argsort(idxperm)

    # Initialize result matrix
    res = np.empty((res_N, msamples), dtype = np.int32)

    # Sample the initial state for all samples
    res[0, :] = np.random.choice(q, size=msamples, p=p0)

    for site in range(N):
        Js = J[site]
        h = H[site]

         # Compute total fields for all samples
        tot_H = np.tile(h, (msamples, 1)).T
        for i in range(site + 1):
            # Vectorized coupling computation
            tot_H += Js[:, res[i, :], i]
        
        # Apply softmax and sample
        p = softmax(tot_H)
        res[site + 1, :] = np.array([np.random.choice(q, p=p[:, i]) for i in range(msamples)])
    
    # Apply permutation
    return res[back_order, :]


def sample_with_weights(arnet: ArNet, msamples: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Return a generated alignment in the form of a N x msamples matrix and the relative probabilities under the module
    

    Args:
        arnet: ArNet object containing model params
        msamples: Number of samples to generate

    Returns:
        Array of shape (N + 1, msamples) of type `::Matrix{Int}` and probabilities
    """ 

    H, J, p0, idxperm = arnet.H, arnet.J, arnet.p0, arnet.idxperm
    q = p0.size
    N, _ = H.shape
    res_N = N + 1

    back_order = np.argsort(idxperm)

    # Initialize result and weight vector
    W = np.empty(msamples, dtype=np.float64)
    res = np.empty((res_N, msamples), dtype=np.int32)

    # Sample initial state
    res[0, :] = np.random.choice(q, size=msamples, p=p0)
    logw = np.log(p0[res[0, :]])

    for site in range(N):
        Js = J[site]
        h = H[site]

        # Compute total fields for all samples
        tot_H = np.tile(h, (msamples, 1)).T
        for i in range(site + 1):
            tot_H += Js[:, res[i, :], i]
        
        p = softmax(tot_H)
        res[site + 1, :] = np.array([np.random.choice(q, p=p[:, i]) for i in range(msamples)])
        logw += np.log(p0[res[site + 1, :]])
    W = np.exp(logw)

    return W, res[back_order, :]