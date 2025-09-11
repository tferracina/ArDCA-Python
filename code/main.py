import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from typing import Tuple, Union

from classes import *
from utils import *


class ArDCA(nn.Module):
    def __init__(self, L: int, q: int, lambda_h: float = 1e-6, lambda_J: float = 1e-4):
        super().__init__()
        self.L, self.q = L, q
        self.lambda_h = lambda_h
        self.lambda_J = lambda_J

        self.h = nn.Parameter(torch.zeros((L, q)))
        self.J = nn.Parameter(torch.zeros((L, L, q, q)))

        tril = torch.tril(torch.ones(L, L, dtype=torch.bool), diagonal=-1)
        self.register_buffer("lower_mask", tril)

    def logits(self, i: int, history: torch.Tensor) -> torch.Tensor:
        """
        Compute unnormalized logits for site i given history a_<i>.
        history: (M, i) tensor of indices for previous residues
        Returns: (M, q) logits
        """
        if i == 0:
            if isinstance(history, int):
                M = history
            else:
                M = history.size(0)
            return self.h[0].unsqueeze(0).expand(M, self.q)
        # i > 0
        H = F.one_hot(history, num_classes=self.q).to(self.h.dtype)  # (M,i,q)
        Ji = self.J[i, :i]                                          # (i,q,q)
        contrib = torch.einsum('jqk,bjk->bq', Ji, H)                # (M,q)
        return self.h[i].unsqueeze(0) + contrib

    def conditional(self, i: int, history: torch.Tensor) -> torch.Tensor:
        """ P(a_i | a_<i>) return (M, q)"""
        return torch.softmax(self.logits(i, history), dim=-1)

    @torch.no_grad()
    def log_prob(self, seqs: torch.Tensor) -> torch.Tensor:
        """
        Compute log-probability of a batch of sequences.
        seq: (M, L)
        Returns: (M,)
        """
        self.eval()
        M, L = seqs.shape
        lp = torch.zeros(M, device=seqs.device, dtype=self.h.dtype)
        for i in range(L):
            logp = torch.log_softmax(self.logits(i, seqs[:, :i]), dim=-1)
            lp += logp.gather(1, seqs[:, i:i+1]).squeeze(1)
        return lp

    @torch.no_grad()
    def energy(self, seqs: torch.Tensor) -> torch.Tensor:
        return -self.log_prob(seqs)
    
    @torch.no_grad()
    def sample(self, n: int = 1, device: str | None = None, rng: torch.Generator | None = None):
        """Autoregressively sample n sequences"""
        self.eval()
        device = device or self.h.device
        seqs = torch.zeros((n, self.L), dtype=torch.long, device=device)
        for i in range(self.L):
            probs = self.conditional(i, seqs[:, :i])
            seqs[:, i] = torch.multinomial(probs, num_samples=1, generator=rng).squeeze(1)
        return seqs

    @torch.no_grad()
    def init_fcolumn_from_freqs(self, seqs: torch.Tensor, eps: float = 1e-8):
        """Lower-triangular mask"""
        M = seqs.size(0)
        f = torch.bincount(seqs[:, 0], minlength=self.q).float()
        f = (f + eps) / (f.sum() + self.q * eps)
        h0 = torch.log(f)
        h0 = h0 - h0.mean()
        self.h[0].copy_(h0)
    

    # TRAINING
    def forward(self, seqs: torch.Tensor, weights: torch.Tensor):
        """
        Compute reweighted NLL + L2 reg
        """
        device = self.h.device
        M, L = seqs.shape
        assert L == self.L
        assert seqs.dtype == torch.long and seqs.shape == (M, L)
        assert (seqs >= 0).all() and (seqs < self.q).all()

        w = weights.to(device)
        M_eff = w.sum().clamp_min(1e-12)
        w_norm = (w / M_eff)

        nll = seqs.new_zeros((), dtype=torch.float32).to(device)

        # i = 0
        logp0 = torch.log_softmax(self.logits(0, M), dim=-1)         # pass M as a tiny convenience
        nll -= (w_norm * logp0.gather(1, seqs[:, 0:1]).squeeze(1)).sum()

        # i >= 1
        for i in range(1, L):
            logp = torch.log_softmax(self.logits(i, seqs[:, :i]), dim=-1)
            nll -= (w_norm * logp.gather(1, seqs[:, i:i+1]).squeeze(1)).sum()

        # L2 regularization (only lower triangle for J)
        reg_h = 0.5 * self.lambda_h * (self.h ** 2).sum()
        J_lower = self.J.masked_select(self.lower_mask[..., None, None])  \
                      .view(-1, self.q, self.q)                           # (#pairs, q, q)
        reg_J = 0.5 * self.lambda_J * (J_lower ** 2).sum()

        loss = nll + reg_h + reg_J
        return loss, nll.detach(), (reg_h + reg_J).detach()

    
def train_ardca_torch(seqs: torch.Tensor,
                      weights: torch.Tensor,
                      q: int,
                      lambda_h: float = 1e-6,
                      lambda_J: float = 1e-4,
                      max_iters: int = 300,
                      lr: float = 1.0,
                      optimizer: str = "LBFGS",
                      freq_h0: bool = True,
                      device: str = "cpu"):
    M, L = seqs.shape
    seqs = seqs.to(device)
    weights = weights.to(device)

    model = ArDCA(L=L, q=q, lambda_h=lambda_h, lambda_J=lambda_J).to(device)
    if freq_h0:
        model.init_fcolumn_from_freqs(seqs)

    if optimizer == "LBFGS":
        opt = torch.optim.LBFGS(model.parameters(), lr=lr, max_iter=max_iters)
        def closure():
            opt.zero_grad(set_to_none=True)
            loss, _, _ = model(seqs, weights)
            loss.backward()
            return loss.detach()
        print(f"Training with LBFGS optimizer (max_iters={max_iters})...")
        opt.step(closure)
    elif optimizer == "Adam":
        opt = torch.optim.Adam(model.parameters(), lr=1e-2)
        print(f"Training with Adam optimizer...")
        for i in range(max_iters):
            opt.zero_grad(set_to_none=True)
            loss, _, _ = model(seqs, weights)
            loss.backward()
            opt.step()
            if (i + 1) % 50 == 0 or i == 0:  # Print every 50 iterations
                print(f"Iteration {i+1}/{max_iters}, Loss: {loss.item():.4f}")
    else:
        raise ValueError("Optimizer not configured. Pick LBFGS or Adam.")
    
    with torch.no_grad():
        loss, nll, reg = model(seqs, weights)
    return model, {"loss": float(loss.cpu()), "nll": float(nll.cpu()), "reg": float(reg.cpu())}


def train_from_fasta(path, max_gap_fraction=1.0, theta=0.8, 
                     lambda_h=1e-6, lambda_J=1e-4, optimizer="LBFGS",
                     device="cpu"):
    print("Reading fasta MSA")
    X = read_fasta_alignment(path, max_gap_fraction)   # (M,L) int8
    print("Computing weights")
    W, M_eff = compute_weights(X, theta, gap_idx=0, count_gaps_as_match=False)
    print("Initializing seqs and weights")
    seqs = torch.from_numpy(X.astype(np.int64)).to(device)
    weights = torch.from_numpy(W.astype(np.float32)).to(device)
    print("Beginning training")
    model, metrics = train_ardca_torch(seqs, weights, q=21,
                                       lambda_h=lambda_h, lambda_J=lambda_J,
                                       optimizer=optimizer, device=device)
    return model, metrics, {"M": X.shape[0], "L": X.shape[1], "M_eff": M_eff}
