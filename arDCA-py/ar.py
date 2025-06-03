import numpy as np
import nlopt
import time
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

def minimize_arnet(alg: ArAlg, var: ArVar) -> Tuple[np.ndarray, np.ndarray]:
    """Site-by-site nonlinear optimization to fit parameters of autoregressive model"""
    N, q, q2, idxperm = var.N, var.q, var.q2, var.idxperm
    epsconv, maxit, method = alg.epsconv, alg.maxit, alg.method
    
    vec_size = (N * (N - 1) // 2) * q2 + (N - 1) * q

    vecps = np.zeros(N - 1, dtype=np.float64) # array of minimal negative log-pseudo likelihoods, one entry per site
    theta = np.zeros(vec_size, dtype=np.float64) # concatenated vector fo all site-wise parameters
    
    for site in range(1, N):
        x0 = np.zeros(site * q2 + q, dtype=np.float64) #site*q2 coupling weights + q local fields

        # construct NLopt optimizer
        opt = nlopt.opt(method, x0.size)

        # set tolerances on absolute/relative function and parameter values
        opt.set_ftol_abs(epsconv)
        opt.set_ftol_rel(epsconv)
        opt.set_xtol_abs(epsconv)
        opt.set_xtol_rel(epsconv)

        # stop after at most maxit evaluations
        opt.set_maxeval(maxit)

        # define objective function
        opt.set_min_objective() # (x, g) -> optimfunwrapper(x, g, site, var) should be minimized
        
        # run and time optimization
        start_time = time.perf_counter()
        minx = opt.optimize(x0)
        elapsed = time.perf_counter() - start_time

        minf = opt.last_optimum_value() # minimized pseudo-likelihood
        status = opt.last_optimize_result()

        if alg.verbose:
            print(f"site = {idxperm[site]}\tpl = {minf:.4f}\ttime = {elapsed:.4f}\t")
            print(f"status: {status}")

        vecps[site - 1] = minf
        offset = (site * (site - 1) // 2) * q2 + (site - 1) * q
        theta[offset : offset + site * q2 + q ] = minx
    
    return theta, vecps


