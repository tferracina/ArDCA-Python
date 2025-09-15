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

    def compute_ar_logits(self, X_oh: torch.Tensor):
        """
        z[m,i,a] = h[i,a] + sum_{j<i} sum_b J[i,j,a,b] * X[m,j,b]
        X_idx: (M, L)
        returns: (M, L, q)
        """
        X_oh = X_oh.to(self.J.dtype)
        M, L, q = X_oh.shape
        logits = self.h_pos.unsqueeze(0).expand(M, -1, -1).clone()  # (M,L,q)

        # Block-sparse: collect lower-triangular indices
        i_idx, j_idx = torch.tril_indices(L, L, offset=-1)
        J_blocks = self.J[i_idx, j_idx]   # (n_pairs, q, q)
        X_blocks = X_oh[:, j_idx]         # (M, n_pairs, q)

        contrib = torch.einsum("mpq,pqr->mpr", X_blocks, J_blocks)  # (M,n_pairs,q)
        logits = logits.index_add(1, i_idx, contrib)

        return logits


    
    def loss(self, X_oh: torch.Tensor, X_idx: torch.Tensor, W: torch.Tensor, lambda_J: Optional[float] = None,
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
        logits = self.compute_ar_logits(X_oh)        # (M,L,q)
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
    def evaluate(self, X_oh: torch.Tensor, X_idx: torch.Tensor, W: torch.Tensor) -> Dict[str, float]:
        """Evaluate model performance on data."""
        self.eval()
        
        logits = self.compute_ar_logits(X_oh)
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
    
    @torch.no_grad()
    def sample(self, n_samples: int = 1, device: str = "cpu") -> np.ndarray:
        """
        Generate new sequences from the trained arDCA model.
        Returns: (n_samples, L) numpy array of sampled sequence indices.
        """
        self.eval()
        L, q = self.L, self.q
        samples = torch.zeros((n_samples, L), dtype=torch.long, device=device)
        X_oh = torch.zeros((n_samples, L, q), dtype=self.J.dtype, device=device)

        for i in range(L):
            logits = self.h_pos[i].clone()
            if i > 0:
                # Add coupling contributions from previous positions
                for j in range(i):
                    prev_a = samples[:, j]
                    J_ij = self.J[i, j]  # shape: (q, q)
                    # Gather the relevant couplings for each sample
                    contrib = J_ij[:, prev_a].T  # shape: (n_samples, q)
                    logits = logits + contrib
            probs = torch.softmax(logits, dim=-1)  # shape: (n_samples, q)
            a_i = torch.multinomial(probs, num_samples=1).squeeze(-1)
            samples[:, i] = a_i
            X_oh[torch.arange(n_samples), i, a_i] = 1.0

        return samples.cpu().numpy()


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
    
    # One-hot encode
    train_X_oh = encode_sequence(train_seqs.cpu().numpy(), q=msa_data.q, device=device)
    val_X_oh = None
    if val_seqs is not None:
        val_X_oh = encode_sequence(val_seqs.cpu().numpy(), q=msa_data.q, device=device)

    model = model.to(device)
    model.init_parameters(train_seqs)
    
    # Optimizer
    if model_params.optimizer.lower() == "lbfgs":
        optimizer = torch.optim.LBFGS(model.parameters(), lr=1e-3)
    elif model_params.optimizer.lower() == "adam":
        optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=0)
    else: raise ValueError(f"Selected optimizer not configured: {model_params.optimizer}");
    
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
    
    for epoch in range(model_params.max_iters):
        model.train()
        optimizer.zero_grad()
            
        loss, info = model.loss(
                train_X_oh, 
                train_seqs,
                train_weights,
                lambda_J=model_params.lambda_J,
                lambda_h=model_params.lambda_h
        )
        loss.backward()
        optimizer.step()
        model.clamp_unused()
        
        # Record training metrics
        history['train_loss'].append(loss.item())
        history['train_nll'].append(info['nll'].item())
        
        # Validation
        if epoch % eval_interval == 0:
            if val_seqs is not None:
                model.eval()
                with torch.no_grad():
                    val_loss, val_info = model.loss(
                        val_X_oh, 
                        val_seqs,
                        val_weights,
                        lambda_J=model_params.lambda_J,
                        lambda_h=model_params.lambda_h
                    )
                    val_metrics = model.evaluate(val_X_oh, val_seqs, val_weights)
                
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
                          optimizer: str = "adam",
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
        optimizer=optimizer,
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

    train_X_oh = encode_sequence(train_seqs, q=msa_data.q, device=device)
    val_X_oh = encode_sequence(val_seqs, q=msa_data.q, device=device)
    train_data = (train_X_oh, torch.tensor(train_seqs, dtype=torch.long).to(device), torch.tensor(train_weights, dtype=torch.float32).to(device))
    val_data = (val_X_oh, torch.tensor(val_seqs, dtype=torch.long).to(device), torch.tensor(val_weights, dtype=torch.float32).to(device))
    
    
    train_metrics = model.evaluate(*train_data)
    val_metrics = model.evaluate(*val_data)
    
    print(f"Final train NLL: {train_metrics['nll_per_pos']:.4f}")
    print(f"Final val NLL: {val_metrics['nll_per_pos']:.4f}")
    print(f"Final val perplexity: {val_metrics['perplexity']:.4f}")

    save_ardca_model(model, "ardca_model.pt")
    
    return model, history, msa_data


def save_ardca_model(model, save_path, extra_params=None):
    """
    Save arDCA model state and parameters.
    """
    checkpoint = {
        'state_dict': model.state_dict(),
        'L': model.L,
        'q': model.q,
        'lambda_h': model.lambda_h,
        'lambda_J': model.lambda_J,
    }
    if extra_params is not None:
        checkpoint.update(extra_params)
    torch.save(checkpoint, save_path)

def load_ardca_model(load_path, device='cpu'):
    """
    Load arDCA model from checkpoint.
    """
    checkpoint = torch.load(load_path, map_location=device)
    model = ArDCA(
        L=checkpoint['L'],
        q=checkpoint['q'],
        lambda_h=checkpoint.get('lambda_h', 1e-6),
        lambda_J=checkpoint.get('lambda_J', 1e-4)
    )
    model.load_state_dict(checkpoint['state_dict'])
    model = model.to(device)
    model.eval()
    return model