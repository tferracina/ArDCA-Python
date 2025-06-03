import numpy as np
from typing import Union, Tuple

from types import ArVar, ArAlg, ArNet
from utils import read_fasta, checkpermorder


def ardca_zw(Z: np.ndarray, W: np.ndarray,
            lambdaJ: float = 0.01,
            lambdaH: float = 0.01,
            pc_factor: float = 0,
            epsconv: float = 1.0e-5,
            maxit: int = 1000,
            verbose: bool = True,
            method: str = "LD_LBFGS",
            permorder: Union[str, np.ndarray] = "ENTROPIC"):

    checkpermorder(permorder)
    
    # Check W is normalized
    if (W < 0).any():
        raise ValueError("weights should be positive")    
    if not np.sum(W) == 1:
        raise ValueError("weights should be normalized")

    Z_copy = np.copy(Z)
    N, M = Z_copy.shape
    if M != len(W):
        raise ValueError("columns in Z should match length of W")

    q = int(Z_copy.max())
    
    #Initialize algorithm parameters in ArAlg object
    aralg = ArAlg(method, verbose, epsconv, maxit)

    #Initialize model variables in ArVar object
    arvar = ArVar(N, M, q, lambdaJ, lambdaH, Z_copy, W, pc_factor, permorder)

    #Perform optimization, calling minimize_arnet
    theta, psval = minimize_arnet(aralg, arvar)

    #Return the result as a tuple of trained model and its parameters
    return ArNet(theta,arvar), arvar


def ardca_fasta(filename: str, **kwargs):
    W, Z, _, _, _ = read_fasta(filename, kwargs.get("max_gap_fraction", 0.9),
                               kwargs.get("theta", "auto"),
                               kwargs.get("remove_dups", True))
    W /= W.sum() # Normalize weights
    return ardca_zw(Z, W, **kwargs)

def minimize_arnet(alg: ArAlg, var: ArVar) -> Tuple[theta, vecps]:
    N, q, q2, idxperm = var.N, var.q, var.q2, var.idxperm
    epsconv, maxit, method = alg.epsconv, alg.maxit, alg.method

    x0 = np.zeros()
