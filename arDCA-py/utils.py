import numpy as np
from Bio import AlignIO

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
    seq_length = alignment.get_alignment_length()
    num_seqs = len(filtered_seqs)
    matrix = np.zeros((seq_length, num_seqs), dtype = np.int8)

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


def read_fasta(filename: str, max_gap_fraction: float, remove_dups: bool): #add theta: Any
    """Read FASTA file"""

    Z = read_fasta_alignment(filename, max_gap_fraction)
    if remove_dups:
        Z = remove_duplicate_sequences(Z)
