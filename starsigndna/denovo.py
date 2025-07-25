"""
Module for de novo mutational signature discovery.

This module provides functionality to discover mutational signatures from observed
mutational data using an expectation-maximization (EM) approach combined with
gradient descent optimization.
"""

from typing import Tuple, Optional
import numpy as np
import warnings
from scipy.stats import poisson
from .main_fixed_denovo import (
    running_simulation_denovo,
    Frobinous,
    convergence,
    cos_sim_matrix
)

def denovo(
    M: np.ndarray,
    n_signatures: int,
    lambd: float,
    O: Optional[np.ndarray] = None,
    em_steps: int = 100,
    gd_steps: int = 50
) -> Tuple[np.ndarray, np.ndarray]:
    """Discover mutational signatures using de novo analysis.
    
    This function performs de novo signature discovery using an expectation-maximization
    approach combined with gradient descent optimization. It attempts to find both the
    signature matrix S and exposure matrix E that best explain the observed mutations M.
    
    Args:
        M: Matrix of observed mutational counts, shape (n_samples, n_mutations)
        n_signatures: Number of signatures to discover
        lambd: Regularization parameter controlling trade-off between fit and sparsity.
               Higher values encourage sparser solutions.
        O: Optional weight matrix indicating mutation potential,
           shape (n_samples, n_mutations) or (n_mutations,)
        em_steps: Maximum number of EM iterations to perform
        gd_steps: Number of gradient descent steps per EM iteration
    
    Returns:
        Tuple[np.ndarray, np.ndarray]: 
            - Exposure matrix E showing signature contributions for each sample,
              shape (n_samples, n_signatures)
            - Signature matrix S containing discovered signatures,
              shape (n_signatures, n_mutations)
    
    Raises:
        Exception: If degenerate signatures are found (rank(S) < len(S))
    
    Notes:
        - The algorithm alternates between optimizing E and S matrices
        - Convergence is checked using the Poisson log-likelihood
        - Both E and S matrices are normalized so rows sum to 1
        - The algorithm may terminate early if convergence is detected
    """
    # Suppress warnings during optimization
    warnings.filterwarnings('ignore')
    
    # Convert inputs to numpy arrays and get dimensions
    M = np.asanyarray(M)
    n_samples, n_mutations = M.shape
    
    # Initialize or validate weight matrix
    if O is None:
        O = np.ones((n_samples, n_mutations), dtype=float)
    M, O = (np.asarray(a) for a in (M, O))
    
    # Initialize exposure matrix with small random values
    E = np.abs(np.random.laplace(loc=0, scale=1, size=(n_samples, n_signatures)))
    E = np.full_like(E, 1e-5)  # Small positive values to avoid numerical issues
    
    # Initialize signature matrix with random values
    S = np.random.rand(n_signatures, n_mutations)
    
    # Optimization parameters
    topt = np.float64("Inf")  
    tedge = np.float64("Inf") 
    
    # Initialize convergence tracking
    pmf_old = np.inf  # Previous log-likelihood
    conv_iter_1 = -1  # First convergence iteration
    conv_check = 0    # Convergence counter
    
    # Main EM loop
    for i in range(em_steps):
        print(f"EM step: {i}")
        
        # Ensure non-negative exposures
        E = np.maximum(E, 0)
        
        # Update exposure matrix E
        E = running_simulation_denovo(E, M, S, O, topt, tedge, lambd, n_steps=gd_steps)
        
        # Update signature matrix S
        S = running_simulation_denovo(S.T, M.T, E.T, O.T, topt, tedge, 0, n_steps=gd_steps).T
        
        # Check for degenerate signatures
        if np.linalg.matrix_rank(S) < S.shape[0]:
            raise Exception(
                "Degenerate signatures found (rank(S) < len(S)). "
                "Try reducing the number of signatures."
            )
        
        # Calculate fit metrics
        mse = Frobinous(M, S, E, O)
        loss = -poisson.logpmf(M, (E @ S) * O)
        mean_loss = np.mean(loss)
        
        # Check convergence
        has_converged = convergence(pmf_old, mean_loss)
        if has_converged:
            print("Convergence detected")
            if conv_iter_1 == -1:
                conv_iter_1 = i
                conv_check = 0
            else:
                conv_check += 1
        else:
            conv_iter_1 = -1
            conv_check = 0
            
        # Early stopping if convergence is stable
        if conv_check == 1:
            print("Algorithm has converged - stopping early")
            break
            
        pmf_old = mean_loss
        print(f"Current negative log-likelihood: {mean_loss:.3e}")
    
    # Post-process results
    E = np.maximum(E, 0)  # Ensure non-negative values
    E /= E.sum(axis=-1, keepdims=True)  # Normalize exposures
    E[np.isnan(E)] = 0  # Handle numerical instabilities
    
    S = np.maximum(S, 0)  # Ensure non-negative values
    S /= S.sum(axis=-1, keepdims=True)  # Normalize signatures
    S[np.isnan(S)] = 0  # Handle numerical instabilities
    
    # Final diagnostics
    final_mse = Frobinous(M, S, E, O)
    print(f"Final MSE: {final_mse:.3e}")
    
    return np.array(E), np.array(S)
