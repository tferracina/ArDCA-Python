import torch
import numpy as np
from dataclasses import dataclass
from typing import Optional, List, Dict


@dataclass
class MSAData:
    seqs: np.ndarray         # [M, L]
    weights: np.ndarray      # [M,], reweighting coeffs w_m
    M_eff: float             # sum(weights)
    L: int
    q: int
    identity_tresh: float


@dataclass
class ModelParams:
    lambda_h: float
    lambda_J: float
    optimizer: str 
    max_iters: int 
    seed: int
    val_frac: float


@dataclass
class TrainState:
    file_path: str
    save_dir: str
    pf: str
    version: int
    lambda_h: float = 1e-6
    lambda_J: float = 1e-4
    max_gap_fraction: float = 0.5
    max_col_gap_fraction: Optional[float] = 0.3
    identity_thresh: float = 0.2
    val_frac: float = 0.1
    max_iters: int = 200
    optimizer: str = 'adam'
    seed: int = 42
    device: str = 'cpu'
