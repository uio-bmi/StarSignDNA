#!/usr/bin/env python3
"""
Training script for mutational signature analysis.

This script performs de novo signature extraction from mutation data using
an alternating optimization approach with regularization.
"""

import numpy as np
from scipy.stats import poisson
import main_fixed_denovo

# Set random seed for reproducibility
np.random.seed(200)


def main(M, n_signatures, lambd):
    """
    Extract mutational signatures from mutation matrix.
    
    Args:
        M (np.ndarray): Mutation matrix (n_samples x n_mutations)
        n_signatures (int): Number of signatures to extract
        lambd (float): Regularization parameter for sparsity
    
    Returns:
        tuple: (exposures, signatures)
            - exposures: Estimated exposure matrix (n_samples x n_signatures)
            - signatures: Estimated signature matrix (n_signatures x n_mutations)
    """
    # Get dimensions
    n_samples, n_mutations = M.shape
    
    # Initialize opportunity matrix (all ones for uniform opportunity)
    O = np.ones((n_samples, n_mutations), dtype=int)
    
    # Initialize exposures with small positive values
    E = np.full((n_samples, n_signatures), 0.00001)
    
    # Initialize signatures with random values
    S = np.random.rand(n_signatures * n_mutations).reshape(n_signatures, n_mutations)
    
    # Set optimization parameters
    topt = np.float64("Inf")
    tedge = np.float64("Inf")
    
    # Ensure non-negative exposures
    if np.any(E < 0):
        E = np.maximum(E, 0)
    
    # Lists to track convergence
    pmf_e = []  # Poisson log-likelihood for exposures
    pmf_s = []  # Poisson log-likelihood for signatures
    d_mse_e = []  # MSE for exposures
    d_mse_s = []  # MSE for signatures
    
    # Initialize convergence tracking
    mse_old = np.inf
    pmf_old = np.inf
    
    # Alternating optimization loop
    for _ in range(25):
        # Update exposures using gradient descent with regularization
        E = main_fixed_denovo.running_simulation_denovo(E, M, S, O, topt, tedge, lambd)
        
        # Update signatures using gradient descent (no regularization)
        S = main_fixed_denovo.running_simulation_denovo(S.T, M.T, E.T, O.T, topt, tedge, 0).T
        
        # Compute MSE for convergence tracking
        mse_e = main_fixed_denovo.Frobinous(M, S, E, O)
        d_mse_s.append(mse_e)
        
        # Compute negative log-likelihood loss
        loss = -poisson.logpmf(M, (E @ S) * O)
        pmf_s.append(np.mean(loss))
        
        # Update convergence tracking
        mse_old = mse_e
        pmf_old = np.mean(loss)
    
    # Post-process exposures: ensure non-negativity and normalize
    if np.any(E < 0):
        E = np.maximum(E, 0)
    E /= E.sum(axis=-1, keepdims=True)
    E[np.isnan(E)] = 0  # Handle NaN values
    
    # Post-process signatures: ensure non-negativity and normalize
    if np.any(S < 0):
        S = np.maximum(S, 0)
    S /= S.sum(axis=-1, keepdims=True)
    S[np.isnan(S)] = 0  # Handle NaN values
    
    return np.array(E), np.array(S)
