# ArDCA-Python


## Abstract

Proteins are fundamental to nearly all biological processes, and understanding the relationship between their amino acid sequences, structures, and functions is a key challenge in computational biology. Direct Coupling Analysis (DCA) has emerged as a powerful statistical approach for inferring residue–residue interactions from the vast archives of protein data, enabling contact prediction and generative modeling of protein families. In this work, we present a Python reimplementation of arDCA, an autoregressive DCA model that factorizes the joint sequence probability distribution into conditionals, to allow exact likelihood maximization and efficient parameter inference without requiring Monte Carlo sampling. Our implementation introduces a block-sparse formulation for computing autoregressive logits, significantly reducing computational overhead and memory usage. We evaluate the model on protein families from Pfam, demonstrating smooth convergence of negative log-likelihood (NLL) and high structural plausibility of generated sequences as measured by AlphaFold’s pLDDT scores. While our implementation remains slower than the original Julia version, it provides an easily expandable and practical framework for protein sequence generation. 

## Structure of Repository

### Written Report
- `LF3199684.typ`:  typst source file
- `LF3199684.pdf`:  compiled pdf file
- `references.bib`:  bibliography for typst

### Code
In the `code` directory:
- `out/`: contains all the produced images and plots
- `ardca.py`:  model, training and evaluation of arDCA
- `classes.py`: dataclass definitions
- `utils.py`: utilities for preprocessing, training, and analyzing results
- `training.ipynb`: training notebook for the model and experiments
- `analyzing_results.ipynb`: notebook used to analyze the model-generated samples
- `requirements.txt`: all the python dependencies

In addition, [ColabFold](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb) was used for images, plots, and pLDDT values.
