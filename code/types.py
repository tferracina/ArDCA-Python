import torch
import numpy as np
from dataclasses import dataclass
from typing import List, Dict

@dataclass
class Alphabet:
    tokens: List[str]
    to_idx: Dict[str, int]
    q: int

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

@dataclass 
class TrainState:
    theta: np.ndarray     # flattened param vector
    value: float          # current objective
    grad_norm: float
    iters: int

@dataclass 
class Metrics: 
    nll: float      # avg reweighted neg log-lh
    reg: float      # L2 penalty value
    total: float    # nll+reg
