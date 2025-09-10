import torch
import numpy as np
from dataclasses import dataclass

#alphabet

#MSA

@dataclass
class ReweightingConfig:
    seqid_threshold: float
    count_gaps_as_match: bool

#ArConfig
@dataclass
class ArConfig:
    L: int
    init_h: str #"zeros" / "emp_logfreq"
    init_J_scale: float
    device: str

@dataclass
class RegularizationConfig:
    lambda_h: float 
    lambda_J: float
    mask_J_in_reg: bool

@dataclass
class OptimConfig:
    name: str
    lr: float
    
