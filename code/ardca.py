import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from typing import Optional, Tuple, Union
from tqdm import tqdm

from classes import *
from utils import *


class ArDCA(nn.Module):
    def __init__(self, L: int, q: int, lambda_h: float = 1e-6, lambda_J: float = 1e-4):
        super().__init__()
        self.L, self.q = L, q
        self.lambda_h = lambda_h
        self.lambda_J = lambda_J

        self.h_pos = nn.Parameter(torch.zeros((L, q)))
        self.J = nn.Parameter(torch.zeros((L, L, q, q)))

        tril = torch.tril(torch.ones(L, L, dtype=torch.bool), diagonal=-1)
        J_mask = tril.reshape((L, L, 1, 1))
        self.register_buffer("J_mask", J_mask)
    
    @torch.no_grad()
    def clamp_unused(self):
        """zero out disallowed J entries for safety."""
        self.J.data *= self.J_mask

    def compute_ar_logits(self, X_idx: torch.Tensor):
        """
        z[m,i,a] = h[i,a] + sum_j<i sum_b J[i,j,a,b] * X[m,j,b]
        X_idx: (M, L)          integer states
        returns: (M, L, q)     z[m,i,a]
        """
        M, L = X_idx.shape
        q = self.q

        X_oh = F.one_hot(X_idx, num_classes=q).to(self.J.dtype)          # (M, L, q)

        pair = torch.einsum('ijab,mjb->mia', self.J * self.J_mask, X_oh)        # (M, L, q)
        h = self.h_pos
        return h.unsqueeze(0) + pair
    
    def loss(self, X_idx: torch.Tensor, W: torch.Tensor, lambda_J: Optional[float] = None,
              lambda_h: Optional[float] = None) -> Tuple[torch.Tensor, Dict[str, torch.Tensor]]:
        """
        Compute weighted negative log-likelihood loss with reg
        X_idx: (M, L) ints
        W:     (M,) weights
        """
        if lambda_J is None:
            lambda_J = self.lambda_J
        if lambda_h is None:
            lambda_h = self.lambda_h

        # 1. Logits -> log-probs
        logits = self.compute_ar_logits(X_idx)        # (M,L,q)
        logp = F.log_softmax(logits, dim=-1)          # (M,L,q)

        # 2. Gather log-probabilities of true aa
        M, L = X_idx.shape
        batch_idx = torch.arange(M, device=X_idx.device)
        pos_idx = torch.arange(L, device=X_idx.device)
        gold_logp = logp[batch_idx[:, None], pos_idx[None, :], X_idx]  # (M,L)

        # 3. weighted negative log likelihood
        data_nll = -(W[:, None] * gold_logp).sum()

        # 4. regularization
        regJ = (self.J.pow(2) * self.J_mask).sum()
        regH = self.h_pos.pow(2).sum()
        
        loss = data_nll + lambda_h * regH + lambda_J * regJ
        info = {
            "nll": data_nll.detach(),
            "regJ": regJ.detach(),
            "regH": regH.detach(),
            "total_loss": loss.detach()
        }
        
        return loss, info
    
    @torch.no_grad()
    def init_h0_from_freqs(self, seqs: torch.Tensor, eps: float = 1e-8):
        """initialize h[0] from empirical frequencies at pos 0"""
        if seqs.dim() != 2:
            raise ValueError("seqs should be (M, L) tensor")
        
        M = seqs.size(0)
        f = torch.bincount(seqs[:, 0], minlength=self.q).float()
        f = (f + eps) / (f.sum() + self.q * eps)
        h0 = torch.log(f)
        h0 = h0 - h0.mean()
        self.h_pos[0].copy_(h0)

    @torch.no_grad()
    def init_parameters(self, seqs: torch.Tensor, init_scale: float = 0.01):
        """Initialize parameters with small random values and set h0 from frequencies."""
        # Initialize h0 from data
        self.init_h0_from_freqs(seqs)
        
        # Small random initialization for other positions
        start_idx = 0
        self.h_pos[start_idx:].normal_(0, init_scale)
        self.J.normal_(0, init_scale)
        
        self.clamp_unused()
    
    def compute_effective_sample_size(self, W: torch.Tensor) -> float:
        """Compute effective sample size from sequence weights."""
        return 1.0 / (W ** 2).sum()
    
    @torch.no_grad()
    def evaluate(self, X_idx: torch.Tensor, W: torch.Tensor) -> Dict[str, float]:
        """Evaluate model performance on data."""
        self.eval()
        
        logits = self.compute_ar_logits(X_idx)
        logp = F.log_softmax(logits, dim=-1)
        
        M, L = X_idx.shape
        batch_idx = torch.arange(M, device=X_idx.device)
        pos_idx = torch.arange(L, device=X_idx.device)
        gold_logp = logp[batch_idx[:, None], pos_idx[None, :], X_idx]
        
        # Weighted metrics
        nll_per_pos = -(W[:, None] * gold_logp)  # (M, L)
        
        metrics = {
            "nll_total": nll_per_pos.sum().item(),
            "nll_per_seq": nll_per_pos.sum(dim=1).mean().item(),
            "nll_per_pos": nll_per_pos.mean().item(),
            "perplexity": torch.exp(nll_per_pos.mean()).item(),
            "effective_sample_size": self.compute_effective_sample_size(W)
        }
        
        return metrics
    
    # @torch.no_grad()
    # def sample(self, n: int = 1, device: str | None = None, rng: torch.Generator | None = None):
    #     """Autoregressively sample n sequences"""
    #     self.eval()
    #     device = device or self.h.device
    #     seqs = torch.zeros((n, self.L), dtype=torch.long, device=device)
    #     for i in range(self.L):
    #         probs = self.conditional(i, seqs[:, :i])
    #         seqs[:, i] = torch.multinomial(probs, num_samples=1, generator=rng).squeeze(1)
    #     return seqs


def train_ardca(model: ArDCA, 
                msa_data: MSAData, 
                model_params: ModelParams,
                device: str = 'cpu') -> Tuple[ArDCA, Dict[str, list]]:
    """
    Train ArDCA model with early stopping.
    model: ArDCA model instance
    msa_data: sequences and weights
    """
    train_seqs, val_seqs = None, None
    train_weights, val_weights = None, None
    
    if model_params.val_frac > 0:
        train_seqs, train_weights, val_seqs, val_weights = split_sequences(
            msa_data.seqs, msa_data.weights, model_params.val_frac, model_params.seed
        )
        # Convert to torch tensors
        train_seqs = torch.tensor(train_seqs, dtype=torch.long).to(device)
        train_weights = torch.tensor(train_weights, dtype=torch.float32).to(device)
        val_seqs = torch.tensor(val_seqs, dtype=torch.long).to(device)
        val_weights = torch.tensor(val_weights, dtype=torch.float32).to(device)
    else:
        train_seqs = torch.tensor(msa_data.seqs, dtype=torch.long).to(device)
        train_weights = torch.tensor(msa_data.weights, dtype=torch.float32).to(device)
    
    model = model.to(device)
    
    model.init_parameters(train_seqs)
    
    # Optimizer
    if model_params.optimizer.lower() == "lbfgs":
        optimizer = torch.optim.LBFGS(model.parameters(), lr=1e-3)
    elif model_params.optimizer.lower() == "adam":
        optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=0)
    else: raise ValueError(f"Selected optimizer not configured: {model_params.optimizer}");
    
    # Data loader
    batch_size = min(1000, len(train_seqs)) 
    train_loader = create_data_loader(train_seqs, train_weights, batch_size, shuffle=True)
    
    # Training history
    history = {
        'train_loss': [],
        'train_nll': [],
        'val_loss': [],
        'val_nll': [],
        'val_perplexity': []
    }
    
    best_val_loss = float('inf')
    patience_counter = 0
    patience = 20
    eval_interval = 10
    
    for epoch in tqdm(range(model_params.max_iters), desc="Training"):
        model.train()
        epoch_loss = 0.0
        epoch_nll = 0.0
        
        for batch_seqs, batch_weights in train_loader:
            if model_params.optimizer.lower() == "lbfgs":
                # Define the closure function for LBFGS
                def closure():
                    optimizer.zero_grad()
                    loss, _ = model.loss(
                        batch_seqs, 
                        batch_weights,
                        lambda_J=model_params.lambda_J,
                        lambda_h=model_params.lambda_h
                    )
                    loss.backward()
                    return loss
                
                # Step with closure for LBFGS
                optimizer.step(closure)
                
                # Compute loss for logging (after optimization)
                loss, info = model.loss(
                    batch_seqs, 
                    batch_weights,
                    lambda_J=model_params.lambda_J,
                    lambda_h=model_params.lambda_h
                )
            else:
                # Standard optimization for Adam and others
                optimizer.zero_grad()
                
                loss, info = model.loss(
                    batch_seqs, 
                    batch_weights,
                    lambda_J=model_params.lambda_J,
                    lambda_h=model_params.lambda_h
                )
                loss.backward()
                
                optimizer.step()

            model.clamp_unused()

            epoch_loss += loss.item()
            epoch_nll += info['nll'].item()
        
        # Record training metrics
        history['train_loss'].append(epoch_loss / len(train_loader))
        history['train_nll'].append(epoch_nll / len(train_loader))
        
        # Validation
        if epoch % eval_interval == 0:
            if val_seqs is not None:
                model.eval()
                with torch.no_grad():
                    val_loss, val_info = model.loss(
                        val_seqs, 
                        val_weights,
                        lambda_J=model_params.lambda_J,
                        lambda_h=model_params.lambda_h
                    )
                    val_metrics = model.evaluate(val_seqs, val_weights)
                
                history['val_loss'].append(val_loss.item())
                history['val_nll'].append(val_info['nll'].item())
                history['val_perplexity'].append(val_metrics['perplexity'])
                
                print(f"Epoch {epoch}: Train Loss={history['train_loss'][-1]:.4f}, "
                      f"Val Loss={val_loss.item():.4f}, Val Perplexity={val_metrics['perplexity']:.4f}")
                
                # Early stopping
                if val_loss.item() < best_val_loss - 1e-6:
                    best_val_loss = val_loss.item()
                    patience_counter = 0
                else:
                    patience_counter += 1
                
                if patience_counter >= patience:
                    print(f"Early stopping at epoch {epoch}")
                    break
            else:
                print(f"Epoch {epoch}: Train Loss={history['train_loss'][-1]:.4f}")
    
    return model, history


def main_training_pipeline(fasta_file: str,
                          lambda_h: float = 1e-6,
                          lambda_J: float = 1e-4,
                          max_gap_fraction: float = 0.5,
                          max_col_gap_fraction: Optional[float] = 0.3,
                          identity_thresh: float = 0.2,
                          val_frac: float = 0.1,
                          max_iters: int = 200,
                          seed: int = 42,
                          device: str = 'cpu') -> Tuple[ArDCA, Dict[str, list], MSAData]:
    """
    Complete training pipeline for ArDCA.
    
    Args:
        fasta_file: Path to MSA FASTA file
        lambda_h: L2 regularization for fields
        lambda_J: L2 regularization for couplings  
        max_gap_fraction: Max fraction of gaps per sequence
        max_col_gap_fraction: Max fraction of gaps per column
        identity_thresh: Sequence identity threshold for reweighting
        val_frac: Fraction of data for validation
        max_iters: Maximum training iterations
        seed: Random seed
        device: 'cuda' or 'cpu'
    
    Returns:
        (trained_model, training_history, msa_data)
    """
    print("Loading MSA data...")
    X_idx = read_fasta_alignment(
        fasta_file, 
        max_gap_fraction=max_gap_fraction,
        max_col_gap_fraction=max_col_gap_fraction,
    )

    weights, M_eff = compute_weights(X_idx=X_idx, theta=identity_thresh, gap_idx=0)

    msa_data = MSAData(
        seqs=X_idx,
        weights=weights,
        M_eff=M_eff,
        L=X_idx.shape[1],
        q=21,
        identity_tresh=identity_thresh
    )
    
    print(f"MSA shape: {msa_data.seqs.shape}")
    print(f"Effective sequences: {msa_data.M_eff:.1f}")
    print(f"Sequence length: {msa_data.L}")
    print(f"Alphabet size: {msa_data.q}")
    
    # Create model parameters
    model_params = ModelParams(
        lambda_h=lambda_h,
        lambda_J=lambda_J,
        optimizer='LBFGS',
        max_iters=max_iters,
        seed=seed,
        val_frac=val_frac
    )

    # Create model
    model = ArDCA(L=msa_data.L, q=msa_data.q, lambda_h=lambda_h, lambda_J=lambda_J)
    
    print("Training ArDCA model...")
    model, history = train_ardca(
        model=model,
        msa_data=msa_data,
        model_params=model_params,
        device=device
    )
    
    # Evaluate final model
    print("Evaluating model...")
    train_seqs, train_weights, val_seqs, val_weights = split_sequences(
        msa_data.seqs, msa_data.weights, val_frac, seed
    )

    train_data = (torch.tensor(train_seqs, dtype=torch.long).to(device),
                 torch.tensor(train_weights, dtype=torch.float32).to(device))
    
    val_data = (torch.tensor(val_seqs, dtype=torch.long).to(device),
               torch.tensor(val_weights, dtype=torch.float32).to(device))
    
    train_metrics = model.evaluate(*train_data)
    val_metrics = model.evaluate(*val_data)
    
    print(f"Final train NLL: {train_metrics['nll_per_pos']:.4f}")
    print(f"Final val NLL: {val_metrics['nll_per_pos']:.4f}")
    print(f"Final val perplexity: {val_metrics['perplexity']:.4f}")
    
    return model, history, msa_data