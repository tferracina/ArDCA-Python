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
        opt.set_min_objective(lambda x, g: compute_pslikeandgrad(x, g, site, var))
        
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

def compute_pslikeandgrad(x: np.ndarray, grad: np.ndarray, site: int, var: ArVar) -> np.float64:
    """ Function that computs the pseudo log-likelihood and gradient
    
    Returns: a Float64 value equal to the (penalized) pseudo-log-likelihood at x.
    Also, writes into grad[:] the gradient ∇ of that same objective.
    """
    N, M, q, q2, lambdaJ, lambdaH, Z, W, IdxZ = var.N, var.M, var.q, var.q2, var.lambdaJ, var.lambdaH, var.Z, var.W, var.IdxZ
    LL = x.shape[0]

    # initialize grad array
    grad = np.zeros_like(x)

    # Apply lambdaJ scaling for indices 1 to LL - q
    grad[:LL - q] = 2.0 * lambdaJ * x[:LL - q]

    # Apply lambdaH scaling for indices LL - q + 1 to LL
    grad[LL - q:] = 2.0 * lambdaH * x[LL - q:]
    # initialize variables
    pseudolike = 0.0 # minus–log-likelihood over all M cases

    vecenergies = np.zeros(q, dtype=np.float64) # vector unnormalized log-potentials, for each state at current site
    expvecenesumnorm = np.zeros(q, dtype=np.float64) # normalized probabilities

    for m in range(M):
        izm = IdxZ[:, m]
        zsm = Z[site, m]

        fill_vecenergies(vecenergies, x, site, izm, q, N)

        lnorm = logsumexp(vecenergies)
        exp_vecenergies_sumnorm = np.exp(vecenergies - lnorm)

        pseudolike -= W[m] * (vecenergies[zsm] - lnorm)

        sq2 = site * q2

        for i in range(site):
            for s in range(q):
                grad[izm[i] + s] += W[m] * exp_vecenergies_sumnorm[s]

            grad[izm[i] + zsm] -= W[m]

        for s in range(q):
            grad[sq2 + s] += W[m] * exp_vecenergies_sumnorm[s]
        
        grad[sq2 + zsm] -= W[m]

    pseudolike += l2norm_asym(x, var)

    return pseudolike


def fill_vecenergies(vecenergies: np.ndarray, x: np.ndarray, site: int, IdxSeq: np.ndarray, q: int, N: int):
    q2 = q * q
    sq2 = site * q2

    # (site × q) array of offsets
    offs = IdxSeq[:, None] + np.arange(q)

    # Sum over axis=0 to get a length-q vector, then add the sq2 term
    vecenergies[:q] = x[offs].sum(axis=0) + x[sq2 : sq2 + q]


def logsumexp(X: np.ndarray):
    u = np.max(X)
    if np.isfinite(u):
        return u + np.log(np.sum(np.exp(X - u)))
    else:
        return u


def l2norm_asym(vec: np.ndarray, var: ArVar):
    q = var.q
    lambdaJ, lambdaH = var.lambdaJ, var.lambdaH

    # slice the coupling-block and the field-block
    J_entries = vec[:-q] # shape ((N-1)*q^2,)
    H_entries = vec[-q:] # shape (q,)

    return lambdaJ * J_entries.dot(J_entries) + lambdaH * H_entries.dot(H_entries)