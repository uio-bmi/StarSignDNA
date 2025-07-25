#!/usr/bin/env python3
"""
Ten-fold cross-validation script for mutational signature analysis.

This script performs signature extraction on training data using specified
hyperparameters (number of signatures k and regularization lambda).
"""

import pandas as pd
import numpy as np
from mutation_sim_train import main

# Extract hyperparameters from snakemake wildcards
k = int(snakemake.wildcards.k)  # Number of signatures
Lambda = float(snakemake.wildcards.Lambda)  # Regularization parameter

# Read the training mutation matrix (exclude the first column which is likely sample IDs)
train_m = pd.read_csv(snakemake.input[0], delimiter="\t")
train_m = train_m.iloc[:, 1:]  # Remove first column (sample IDs)

# Convert to numpy array for processing
train_m = train_m.to_numpy().astype(float)

# Run signature extraction to get exposures and signatures
E, S = main(train_m, k, Lambda)

# Save results
np.savetxt(snakemake.output[0], E)  # Save exposures
np.savetxt(snakemake.output[1], S)  # Save signatures
