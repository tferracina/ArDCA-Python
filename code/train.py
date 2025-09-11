from main import *
from utils import *

fasta_path = "../data/PF00014_mgap6.fasta.gz"

train_from_fasta(fasta_path, max_gap_fraction=1.0, theta=0.8, 
                     lambda_h=1e-6, lambda_J=1e-4, optimizer="LBFGS",
                     device="cpu")