# MSA Tensor 
import torch
import numpy as np


X_idx = [] #shape [B, L]
X_onehot = [] # shape [B, L, q], X[b, i, a] = 1 iff a at i