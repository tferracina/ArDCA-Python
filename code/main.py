# MSA Tensor 
import torch
import torch.nn.functional as F
import numpy as np
from typing import Tuple, Union

# define the model
def model(h: np.ndarray, J: np.ndarray):
    # fields h[i, a]
    # couplings J[i, j, a, b]
