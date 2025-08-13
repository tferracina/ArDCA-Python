import numpy as np

from ar_types import ArNet, ArVar

def epistatic_score(arnet: ArNet, arvar: ArVar, seqid: int, pc: np.float64, min_separation: int):
    H, J, p0, idxperm = arnet.H, arnet.J, arnet.p0, arnet.idxperm
    Z, M, N, q = arvar.Z, arvar.M, arvar.N, arvar.q

    if not (1 <= seqid <= M):
        raise ValueError(f"seqid {seqid} should be in the interval [1, ... , {M}]")
    
    Da = np.empty((q, N))
    Dab = np.empty((q, q, N, N))

    xori = Z[:, seqid]

    ppc = (1-pc) * p0 + pc * np.ones(q)/q
