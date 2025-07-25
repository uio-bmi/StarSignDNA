"""
Module for refitting mutational signatures to observed data.

This module provides functionality to refit known mutational signatures to observed
mutational data using a regularized optimization approach.
"""

from typing import Tuple, Optional
import numpy as np
import warnings
import math
from scipy.stats import poisson, entropy
from numpy import count_nonzero
from .main_fixed_denovo import (
    Frobinous,
    Frobinous_reconstuct,
    running_simulation_refit
)
from scipy import stats

def refit(
    M: np.ndarray,
    S: np.ndarray,
    O: Optional[np.ndarray] = None,
    lambd: float = int,
    n_iterations: int = 1000
) -> np.ndarray:
    """Refit known mutational signatures to observed mutational data.
    
    This function performs signature refitting using a regularized optimization approach.
    It attempts to find the optimal exposure matrix E that explains the observed mutations M
    given the signature matrix S, while maintaining sparsity through L1 regularization.
    
    Args:
        M: Matrix of observed mutational counts, shape (n_samples, n_mutations)
        S: Matrix of known mutational signatures, shape (n_signatures, n_mutations)
        O: Optional weight matrix indicating mutation potential,
           shape (n_samples, n_mutations) or (n_mutations,)
        lambd: Regularization parameter controlling trade-off between fit and sparsity.
               Higher values encourage sparser solutions.
        n_iterations: Number of optimization iterations to perform
    
    Returns:
        np.ndarray: Exposure matrix E showing signature contributions for each sample,
                   shape (n_samples, n_signatures)
    
    Notes:
        - The weight matrix O will be tiled to match M's shape if provided as 1D
        - The exposure matrix is normalized so each sample's exposures sum to 1
        - Very small exposure values (<1e-5) are set to 0 for sparsity
        - The function prints diagnostic information about the fit quality
    """
    # Suppress warnings during optimization
    warnings.filterwarnings('ignore')
    np.random.seed(1000)  # For reproducibility
    
    # Get dimensions
    n_samples = len(M)
    n_signatures = len(S)
    n_mutations = 96  # Standard number of mutation types
    
    # Initialize or validate weight matrix
    if O is None:
        O = np.ones((n_samples, n_mutations), dtype=float)
    elif O.shape != M.shape:
        O = np.tile(O, (M.shape[0], 1))
    
    # Convert inputs to numpy arrays
    M, S, O = (np.asarray(a) for a in (M, S, O))
    
    # Initialize exposure matrix with random values
    E = np.abs(np.random.laplace(loc=0, scale=1, size=(n_samples, n_signatures)))
    E = np.full_like(E, 1.0)  # Start with uniform exposures
    
    # Optimization parameters
    topt = np.float64("Inf") 
    tedge = np.float64("Inf") 
    
    # Perform optimization
    E = running_simulation_refit(
        E=E,
        M=M,
        S=S,
        O=O,
        topt=topt,
        tedge=tedge,
        lambd=lambd,
        n_steps=n_iterations
    )
    
    # Post-process results
    E = np.maximum(E, 1e-5)  # Ensure non-negative exposures
    E /= E.sum(axis=-1, keepdims=True)  # Normalize exposures
    E[np.isnan(E)] = 0  # Handle any numerical instabilities
    
    # Calculate fit diagnostics
    mse = Frobinous(M, S, E, O)  # Mean squared error
    loss = -poisson.logpmf(M, (E @ S) * O)  # Negative log likelihood
    loss[loss == np.float64("Inf")] = 0  # Handle infinite values
    sparsity = 1.0 - (count_nonzero(E) / float(E.size))  # Fraction of zero entries
    
    # Print diagnostics
    print(f"Sparsity: {sparsity:.3f}")
    print(f"Mean Squared Error: {np.mean(mse):.3e}")
    print(f"Mean Negative Log Likelihood: {np.mean(loss):.3e}")
    
    return E





