#!/usr/bin/env python3
"""
Train-test split script for mutational signature analysis.

This script splits the input mutation matrix into training and test sets
for cross-validation purposes. It creates 5 different train-test splits
with a 70-30 ratio.
"""

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split

# Read the input mutation matrix
M = pd.read_csv(snakemake.input[0], delimiter="\t")

# Create 5 different train-test splits for cross-validation
for i in range(5):
    # Split the data with 70% training, 30% test
    train, test = train_test_split(M, test_size=0.3, random_state=i)
    
    # Save training and test sets
    train.to_csv(snakemake.output.trains[i], sep="\t", header=True, index=False)
    test.to_csv(snakemake.output.tests[i], sep="\t", header=True, index=False)
