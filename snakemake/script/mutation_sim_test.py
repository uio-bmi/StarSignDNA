#!/usr/bin/env python3
"""
Test script for mutational signature analysis.

This script evaluates the performance of learned signatures by computing
exposures on test data and calculating the negative log-likelihood loss.
"""

import numpy as np
from scipy.stats import poisson
import main_fixed_denovo

# Set random seed for reproducibility
np.random.seed(100)


def main(M, S):
    """
    Evaluate learned signatures on test data.
    
    Args:
        M (np.ndarray): Test mutation matrix (n_samples x n_mutations)
        S (pd.DataFrame): Learned signatures matrix (n_signatures x n_mutations)
    
    Returns:
        tuple: (exposures, mean_loss)
            - exposures: Estimated exposure matrix
            - mean_loss: Mean negative log-likelihood loss
    """
    # Convert signatures to numpy array
    S = S.to_numpy()
    
    # Get dimensions
    n_samples, n_mutations = M.shape
    n_signatures = S.shape[0]
    
    # Initialize opportunity matrix (all ones for uniform opportunity)
    O = np.ones((n_samples, n_mutations), dtype=int)
    
    # Initialize true exposures with Laplace distribution
    E_true = np.abs(np.random.laplace(loc=0, scale=1, 
                                     size=n_samples * n_signatures).reshape(n_samples, n_signatures))
    
    # Initialize estimated exposures with small positive values
    E = np.full_like(E_true, 0.00001)
    
    # Set optimization parameters
    topt = np.float64("Inf")
    tedge = np.float64("Inf")
    
    # Ensure non-negative exposures
    if np.any(E < 0):
        E = np.maximum(E, 0)
    
    # Run signature extraction algorithm
    E = main_fixed_denovo.running_simulation_denovo(E, M, S, O, topt, tedge, 0)
    
    # Compute negative log-likelihood loss using Poisson distribution
    loss = -poisson.logpmf(M, (E @ S) * O)
    loss[loss == np.float64("Inf")] = 0  # Handle infinite values
    
    # Ensure non-negative exposures and normalize
    if np.any(E < 0):
        E = np.maximum(E, 0)
    E /= E.sum(axis=-1, keepdims=True)
    E[np.isnan(E)] = 0  # Handle NaN values
    
    return E, np.mean(loss)
