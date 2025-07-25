#!/usr/bin/env python3
"""
Evaluation script for mutational signature analysis.

This script evaluates the performance of signature extraction by testing
the learned signatures on a test dataset and computing the loss.
"""

import pandas as pd
import numpy as np
from mutation_sim_test import main

# Read the test mutation matrix (exclude the first column which is likely sample IDs)
test_m = pd.read_csv(snakemake.input[0], delimiter="\t")
test_m = test_m.iloc[:, 1:]  # Remove first column (sample IDs)

# Convert to numpy array for processing
test_m = test_m.to_numpy().astype(float)

# Read the learned signatures
S = pd.read_csv(snakemake.input[1], header=None, delimiter=" ")

# Run evaluation to get exposures and loss
E, loss = main(test_m, S)

# Save results
np.savetxt(snakemake.output[0], E)
with open(snakemake.output[1], "w") as f:
    f.write(str(loss))
