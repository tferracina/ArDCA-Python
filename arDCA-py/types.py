from dataclasses import dataclass
from typing import List, Union, Optional
import numpy as np
from numpy.typing import NDArray
import numpy.random as npr

@dataclass
class ArVar:  #holds data and model params
    N: int  # length of sequence
    M: int  # number of sequences
    q: int  # maximum character in sequence
    q2: int
    lambdaJ: float  # Coupling strength parameter
    lambdaH: float  # Field strength parameter
    Z: NDArray  # MSA matrix
    W: NDArray  # Weights vector
    pc: float = 0 # pseudocount factor for p0, defaults to 1/M
    IdxZ: Optional[NDArray] = None  # partial index computation to speed up energy calculation
    idxperm: Optional[NDArray] = None  # permutation index

    # insert ArVar methods


@dataclass
class ArAlg:  # stores algorithm config
    method: str
    verbose: bool
    epsconv: float
    maxit: int

@dataclass
class ArNet:  # implements ar network functionality
    idxperm: NDArray
    p0: NDArray
    J: List[NDArray]  # List of 3D arrays
    H: List[NDArray]  # List of 1D arrays
