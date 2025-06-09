import numpy as np, nlopt, time
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
    if np.any(W <= 0):
        raise ValueError("weights must be positive")    
    if not np.isclose(W.sum(), 1.0, rtol=1e-10):
        raise ValueError("weights should be normalized to 1")

    Z_copy = Z.copy()
    N, M = Z_copy.shape
    if M != W.size:
        raise ValueError("columns in Z should match length of W")

    q = int(Z_copy.max())
    
    #Initialize algorithm parameters in ArAlg object
    aralg = ArAlg(method, verbose, epsconv, maxit)

    #Initialize model variables in ArVar object
    arvar = ArVar(N, M, q, lambdaJ, lambdaH, Z_copy, W, pc_factor, permorder)

    #Perform optimization, calling minimize_arnet
    theta, psval = minimize_arnet(aralg, arvar)
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
    
    vec_size = (N * (N - 1) // 2) * q2 + (N - 1) * q
    vecps = np.empty(N - 1) # array of minimal negative log-pseudo likelihoods, one entry per site
    theta = np.empty(vec_size) # concatenated vector fo all site-wise parameters
    
    for site in range(1, N):
        x0 = np.zeros(site * q2 + q, dtype=np.float64) #site*q2 coupling weights + q local fields

        # construct NLopt optimizer
        opt = nlopt.opt(alg.method, x0.size)

        # set tolerances on absolute/relative function and parameter values
        opt.set_ftol_abs(alg.epsconv); opt.set_ftol_rel(alg.epsconv)
        opt.set_xtol_abs(alg.epsconv); opt.set_xtol_rel(alg.epsconv)

        # stop after at most maxit evaluations
        opt.set_maxeval(alg.maxit)

        # define objective function
        opt.set_min_objective(
            lambda x, g, s=site: compute_pslikeandgrad(x, g, s, var)
            )
        
        # run and time optimization
        start_time = time.process_time()
        minx = opt.optimize(x0)
        elapsed = time.process_time() - start_time
        minf = opt.last_optimum_value() 

        status = opt.last_optimize_result()

        if alg.verbose:
            print(f"site = {idxperm[site]}\tpl = {minf:.4f}\ttime = {elapsed:.4f}\t")
            print(f"status: {status}")

        vecps[site - 1] = minf
        offset = (site * (site - 1) // 2) * q2 + (site - 1) * q
        theta[offset : offset + minx.size ] = minx
    
    return theta, vecps

def compute_pslikeandgrad(x: np.ndarray, g: np.ndarray, site: int, var: ArVar) -> float:
    """ Function that computs the pseudo log-likelihood and gradient
    
    Returns: a Float64 value equal to the (penalized) pseudo-log-likelihood at x.
    Also, writes into grad[:] the gradient ∇ of that same objective.
    """
    # initialize grad array
    if g.size:
        g[:] = 0.0
        grad = g
    else:
        grad = np.zeros_like(x)
    
    N, M, q, q2 = var.N, var.M, var.q, var.q2
    lambdaJ, lambdaH = var.lambdaJ, var.lambdaH
    Z, W, IdxZ = var.Z, var.W, var.IdxZ

    LL = x.shape[0]

    # Apply lambdaJ scaling for indices 1 to LL - q
    grad[:LL - q] = 2.0 * lambdaJ * x[:LL - q]
    # Apply lambdaH scaling for indices LL - q + 1 to LL
    grad[LL - q:] = 2.0 * lambdaH * x[LL - q:]

    # initialize variables
    plike = 0.0 # minus–log-likelihood over all M cases
    sq2 = site * q2
    vecE = np.empty(q) # vector unnormalized log-potentials, for each state at current site

    for m in range(M):
        izm = IdxZ[:site, m]
        zsm = Z[site, m]

        vecE[:] = x[izm].reshape(site, q).sum(0) + x[sq2 : sq2 + q]

        lnorm = logsumexp(vecE)
        probs = np.exp(vecE - lnorm)

        plike -= W[m] * (vecE[zsm] - lnorm)

        for i in range(site):
            grad[izm[i] + izm[i] + q] += W[m] * probs
            grad[izm[i] + zsm]        -= W[m]

        grad[sq2 : sq2 + q] += W[m] * probs
        grad[sq2 + zsm]     -= W[m]
    
    # add L2 penalty
    J = x[: LL - q]
    H = x[LL - q :]
    plike += lambdaJ * J.dot(J) + lambdaH * H.dot(H)

    return plike


def logsumexp(X: np.ndarray):
    u = np.max(X)
    return u + np.log(np.exp(X - u).sum()) if np.isfinite(u) else u