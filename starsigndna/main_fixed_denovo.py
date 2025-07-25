"""
Core optimization functions for de novo mutational signature discovery.

This module provides the core optimization routines used in de novo signature discovery,
including gradient computation, Newton-Raphson optimization, and convergence checking.
"""

from typing import Optional, Tuple, List, Union
import numpy as np
from numpy.random.mtrand import poisson
from scipy.stats import poisson, entropy
import scipy.spatial as sp
import pandas as pd
from scipy.optimize import linear_sum_assignment
import string


def compute_local_gradients(E: np.ndarray, M: np.ndarray, S: np.ndarray, O: np.ndarray) -> np.ndarray:
    """Compute local gradients for refitting exposure optimization.
    
    Args:
        E: Exposure matrix, shape (n_samples, n_signatures)
        M: Mutation count matrix, shape (n_samples, n_mutations)
        S: Signature matrix, shape (n_signatures, n_mutations)
        O: Weight matrix, shape (n_samples, n_mutations)
    
    Returns:
        np.ndarray: Local gradients, shape (n_samples, n_signatures)
    """
    n_samples, n_signatures, n_mutations = (E.shape[0], S.shape[0], M.shape[1])
    local_gradients = np.empty_like(E)
    
    for i in range(n_samples):
        for r in range(n_signatures):
            numerator = M[i] * S[r]
            denominator_sum = np.array([E[i] @ S[:, k] for k in range(n_mutations)])
            denominator_sum_c = denominator_sum + 1e-6  # Numerical stability
            local_gradients[i, r] = np.sum((numerator/denominator_sum_c) - O[i]*S[r])
            
    return local_gradients


def compute_local_gradients_denovo(E: np.ndarray, M: np.ndarray, S: np.ndarray, O: np.ndarray) -> np.ndarray:
    """Compute local gradients for de novo signature discovery.
    
    Similar to compute_local_gradients but without numerical stability term.
    
    Args:
        E: Exposure matrix, shape (n_samples, n_signatures)
        M: Mutation count matrix, shape (n_samples, n_mutations)
        S: Signature matrix, shape (n_signatures, n_mutations)
        O: Weight matrix, shape (n_samples, n_mutations)
    
    Returns:
        np.ndarray: Local gradients, shape (n_samples, n_signatures)
    """
    n_samples, n_signatures, n_mutations = (E.shape[0], S.shape[0], M.shape[1])
    local_gradients = np.empty_like(E)
    
    for i in range(n_samples):
        for r in range(n_signatures):
            numerator = M[i] * S[r]
            denominator_sum = np.array([E[i] @ S[:, k] for k in range(n_mutations)])
            local_gradients[i, r] = np.sum((numerator / denominator_sum) - O[i] * S[r])
            
    return local_gradients


def compute_hessians(E: np.ndarray, M: np.ndarray, S: np.ndarray) -> np.ndarray:
    """Compute Hessian matrices for Newton-Raphson optimization.
    
    Args:
        E: Exposure matrix, shape (n_samples, n_signatures)
        M: Mutation count matrix, shape (n_samples, n_mutations)
        S: Signature matrix, shape (n_signatures, n_mutations)
    
    Returns:
        np.ndarray: Hessian matrices, shape (n_samples, n_signatures, n_signatures)
    """
    denominator = (E @ S) ** 2 + 1e-6  # Numerical stability
    numerator = M[:, None, None, :] * S[None, :, None, :] * S[None, None, :, :]
    res = numerator / denominator[:, None, None, :]
    hessians = -res.sum(axis=-1)
    return hessians


def compute_global_gradient(E: np.ndarray, local_gradients: np.ndarray, lambd: float) -> np.ndarray:
    """Compute global gradient with L1 regularization.
    
    Args:
        E: Exposure matrix, shape (n_samples, n_signatures)
        local_gradients: Local gradients, shape (n_samples, n_signatures)
        lambd: L1 regularization parameter
    
    Returns:
        np.ndarray: Global gradients, shape (n_samples, n_signatures)
    """
    cond_a = local_gradients - lambd * np.sign(E)
    cond_b = local_gradients - lambd * np.sign(local_gradients)
    cond_c = 0
    return np.where(E != 0, cond_a, np.where(np.abs(local_gradients) > lambd, cond_b, cond_c))


def compute_global_gradient_denovo(E: np.ndarray, local_gradients: np.ndarray, 
                                 matrice_lambda: float) -> np.ndarray:
    """Compute global gradient for de novo signature discovery.

    Args:
        E: Exposure matrix, shape (n_samples, n_signatures)
        local_gradients: Local gradients, shape (n_samples, n_signatures)
        matrice_lambda: L1 regularization parameter matrix

    Returns:
        np.ndarray: Global gradients, shape (n_samples, n_signatures)
    """
    cond_a = local_gradients - matrice_lambda * np.sign(E)
    cond_b = local_gradients - matrice_lambda * np.sign(local_gradients)
    cond_c = 0
    return np.where(E != 0, cond_a, np.where(np.all(np.abs(local_gradients) > matrice_lambda), cond_b, cond_c))


def compute_topt(E: np.ndarray, local_gradients: np.ndarray, 
                global_gradients: np.ndarray, hessians: np.ndarray) -> float:
    """Compute optimal step size using Newton-Raphson method.
    
    Args:
        E: Exposure matrix, shape (n_samples, n_signatures)
        local_gradients: Local gradients, shape (n_samples, n_signatures)
        global_gradients: Global gradients, shape (n_samples, n_signatures)
        hessians: Hessian matrices, shape (n_samples, n_signatures, n_signatures)
    
    Returns:
        float: Optimal step size
    """
    if np.any(E < 0):
        E = np.maximum(E, 0)
        
    numerator = np.linalg.norm(global_gradients, ord='fro')
    gg_vectors = (gg[:, None] for gg in global_gradients)
    denominator = sum(gg.T @ h @ gg for gg, h in zip(gg_vectors, hessians))
    
    return -(numerator / denominator) + 1e-5


def compute_t_edge(E: np.ndarray, global_gradients: np.ndarray) -> float:
    """Compute maximum allowable step size to maintain non-negativity.
    
    Args:
        E: Exposure matrix, shape (n_samples, n_signatures)
        global_gradients: Global gradients, shape (n_samples, n_signatures)
    
    Returns:
        float: Maximum allowable step size
    """
    global_gradients = global_gradients.flatten()
    E = E.flatten()
    mask = (np.sign(E) == -np.sign(global_gradients)) & (np.sign(E) != 0)
    
    if not np.any(mask):
        return np.inf
        
    return np.min(-(E/global_gradients)[mask]) + 1e-2


def compute_topt_denovo(E: np.ndarray, local_gradients: np.ndarray,
                       global_gradients: np.ndarray, hessians: np.ndarray) -> float:
    """Compute optimal step size for de novo signature discovery.
    
    Args:
        E: Exposure matrix, shape (n_samples, n_signatures)
        local_gradients: Local gradients, shape (n_samples, n_signatures)
        global_gradients: Global gradients, shape (n_samples, n_signatures)
        hessians: Hessian matrices, shape (n_samples, n_signatures, n_signatures)
    
    Returns:
        float: Optimal step size
    """
    numerator = np.sum(global_gradients * local_gradients)
    gg_vectors = (gg[:, None] for gg in global_gradients)
    denominator = sum(gg.T @ h @ gg for gg, h in zip(gg_vectors, hessians))

    return -(numerator / denominator)


def compute_t_edge_denovo(E: np.ndarray, global_gradients: np.ndarray) -> float:
    """Compute maximum allowable step size for de novo discovery.

    Args:
        E: Exposure matrix, shape (n_samples, n_signatures)
        global_gradients: Global gradients, shape (n_samples, n_signatures)

    Returns:
        float: Maximum allowable step size
    """
    global_gradients = global_gradients.flatten()
    E = E.flatten()
    
    # Remove zero gradient entries
    ind = np.where(global_gradients == 0)
    E_conv = np.delete(E, ind[0])
    global_gradients_conv = np.delete(global_gradients, ind[0])
    
    mask = (np.sign(E_conv) == -np.sign(global_gradients_conv)) & (np.sign(E_conv) != 0)
    
    if not np.any(mask):
        return np.inf
        
    assert np.all(global_gradients_conv != 0)
    return np.min(-(E_conv / global_gradients_conv)[mask])


def min_topt_tedge(topt: float, tedge: float) -> float:
    """Select minimum step size between optimal and edge steps.
    
    Args:
        topt: Optimal step size from Newton-Raphson
        tedge: Maximum allowable step size
    
    Returns:
        float: Minimum of the two step sizes
    """
    topt_tedge = np.minimum(float(topt), float(tedge))
    return topt_tedge


def update_exposure_gradient(E: np.ndarray, global_gradients: np.ndarray, 
                           topt_tedge: float) -> np.ndarray:
    """Update exposures using gradient descent step.
    
    Args:
        E: Exposure matrix, shape (n_samples, n_signatures)
        global_gradients: Global gradients, shape (n_samples, n_signatures)
        topt_tedge: Step size to use
    
    Returns:
        np.ndarray: Updated exposure matrix
    """
    if np.any(E < 0):
       E = np.maximum(E, 0)
    return E + topt_tedge * global_gradients


def newton_raphson1(E: np.ndarray, global_gradients: np.ndarray, 
                   hessians: np.ndarray) -> Optional[np.ndarray]:
    """Perform Newton-Raphson optimization step.
    
    Args:
        E: Exposure matrix, shape (n_samples, n_signatures)
        global_gradients: Global gradients, shape (n_samples, n_signatures)
        hessians: Hessian matrices, shape (n_samples, n_signatures, n_signatures)
    
    Returns:
        Optional[np.ndarray]: Updated exposure matrix, or None if optimization fails
    """
    assert np.all(E >= 0), "Error: E matrix contains negative values"
    nr = []
    
    for E_row, gradient_row, hessian in zip(E, global_gradients, hessians):
        non_active = (E_row == 0) | (gradient_row == 0)
        active_mask = ~non_active
        active_gradients = gradient_row[active_mask]
        active_hessian = hessian[active_mask][:, active_mask]
        new_row = E_row.copy()
        
        det = np.linalg.det(active_hessian)
        if det < 1e-10:
            return None
            
        new_row[active_mask] = E_row[active_mask] - np.linalg.inv(active_hessian) @ active_gradients
        nr.append(new_row)
        
    return np.array(nr)


def update_exposure_NR(E: np.ndarray, global_gradients: np.ndarray,
                      topt: float, tedge: float, new_E: np.ndarray) -> np.ndarray:
    """Update exposures using Newton-Raphson or gradient descent.
    
    Args:
        E: Current exposure matrix
        global_gradients: Global gradients
        topt: Optimal step size
        tedge: Maximum allowable step size
        new_E: Proposed new exposure matrix from Newton-Raphson
    
    Returns:
        np.ndarray: Updated exposure matrix
    """
    if np.any(E < 0):
        E = np.maximum(E, 0)
        assert topt <= tedge
    
    return np.where(np.sign(new_E) == np.sign(E),
                    new_E,
                    E + topt * global_gradients)


def convergence(E: float, E_hat: float, tol: float = 1e-6) -> bool:
    """Check for convergence between iterations.
    
    Args:
        E: Current value
        E_hat: Previous value
        tol: Convergence tolerance
    
    Returns:
        bool: True if converged, False otherwise
    """
    return np.abs((E_hat - E) / E) < tol


def mean_exposure(E: np.ndarray) -> float:
    """Calculate mean exposure across all signatures.
    
    Args:
        E: Exposure matrix
    
    Returns:
        float: Mean exposure value
    """
    return np.mean(E)


def Frobinous(M: np.ndarray, S: np.ndarray, E: np.ndarray, O: np.ndarray) -> float:
    """Compute Frobenius norm of reconstruction error.
    
    Args:
        M: Mutation count matrix
        S: Signature matrix
        E: Exposure matrix
        O: Weight matrix
    
    Returns:
        float: Frobenius norm of difference
    """
    return np.linalg.norm(M - (E @ S) * O, ord='fro')


def Frobinous_reconstuct(M: np.ndarray, S: np.ndarray, E: np.ndarray, 
                        O: np.ndarray) -> Tuple[float, np.ndarray]:
    """Compute Frobenius norm and reconstructed mutation matrix.
    
    Args:
        M: Mutation count matrix
        S: Signature matrix
        E: Exposure matrix
        O: Weight matrix
    
    Returns:
        Tuple[float, np.ndarray]: 
            - Frobenius norm of difference
            - Reconstructed mutation matrix
    """
    M_hat = (E @ S) * O
    M_hat *= M.sum(axis=1).reshape(-1, 1)
    return np.linalg.norm(M - M_hat), M_hat


def sparsity(E: np.ndarray) -> bool:
    """Calculate sparsity of exposure matrix.
    
    Args:
        E: Exposure matrix
    
    Returns:
        bool: True if sparsity >= 0.8, False otherwise
    """
    sparsity = 1.0 - (np.count_nonzero(E) / float(E.size))
    return sparsity >= 0.8


def mse(E: np.ndarray, E_hat: np.ndarray) -> float:
    """Calculate mean squared error between matrices.
    
    Args:
        E: True matrix
        E_hat: Predicted matrix
    
    Returns:
        float: Mean squared error
    """
    return np.square(np.subtract(E, E_hat)).mean()


def cos_sim_matrix(matrix1: np.ndarray, matrix2: pd.DataFrame) -> Tuple[pd.DataFrame, str, pd.DataFrame]:
    """Compute cosine similarity between signature matrices.
    
    Args:
        matrix1: First signature matrix
        matrix2: Second signature matrix with index labels
    
    Returns:
        Tuple[pd.DataFrame, str, pd.DataFrame]:
            - Full cosine similarity matrix
            - String describing optimal assignment
            - Filtered cosine similarity matrix for best matches
    """
    cosine = 1 - sp.distance.cdist(matrix1, matrix2, 'cosine')
    m = 1 - sp.distance.cdist(matrix1, matrix2, 'cosine')
    match_coresp = linear_sum_assignment(m, maximize=True)
    
    match_coresp_1 = "Optimal assignment: \n"
    select_row = []
    for i, j in zip(match_coresp[0], match_coresp[1]):
        match_coresp_1 = match_coresp_1 + str(i) + " => " + str(j) + ", "
        select_row.append(j)
        
    cosine_transpose = np.transpose(cosine)
    cosine_extract = cosine_transpose[select_row]
    
    rows_names = []
    index_names = matrix2.index.values.tolist()
    for index in select_row:
        rows_names.append(index_names[index])
        
    alphabet = list(string.ascii_uppercase)
    cosine_extract_columns = ['Denovo ' + alphabet[k] for k in range(cosine_extract.shape[1])]
    
    df = pd.DataFrame(cosine_extract, columns=cosine_extract_columns)
    df[" "] = rows_names
    df = df.set_index([" "])
    
    cosine_index = pd.DataFrame(cosine_transpose, columns=cosine_extract_columns)
    cosine_index[" "] = index_names
    cosine_index_detect = cosine_index.set_index([" "])
    
    return cosine_index_detect, match_coresp_1[:-2], df


def running_simulation_refit(E: np.ndarray, M: np.ndarray, S: np.ndarray, O: np.ndarray,
                           topt: float, tedge: float, lambd: float, n_steps: int) -> np.ndarray:
    """Run optimization for signature refitting.
    
    Args:
        E: Initial exposure matrix
        M: Mutation count matrix
        S: Signature matrix
        O: Weight matrix
        topt: Initial optimal step size
        tedge: Initial edge step size
        lambd: L1 regularization parameter
        n_steps: Maximum number of optimization steps
    
    Returns:
        np.ndarray: Optimized exposure matrix
    """
    old_loss = np.inf
    loss_hat = np.inf
    conv_iter_1 = -1
    conv_check = 0
    
    for step in range(n_steps):
        print(f"Iteration step: {step}")
        
        if np.any(E < 0):
            E = np.maximum(E, 0)
            
        # Compute gradients and step sizes
        local_gradients = compute_local_gradients(E, M, S, O)
        hessians = compute_hessians(E, M, S)
        global_gradients = compute_global_gradient(E, local_gradients, lambd)
        tedge = compute_t_edge(E, global_gradients)
        topt = compute_topt(E, local_gradients, global_gradients, hessians)
        t_min = min_topt_tedge(topt, tedge)
        
        # Update exposures
        E = update_exposure_gradient(E, global_gradients, t_min)
        
        if topt >= tedge:
            # Gradient descent update
            E = update_exposure_gradient(E, global_gradients, t_min)
        else:
            # Try Newton-Raphson update
            newton_raphason = newton_raphson1(E, global_gradients, hessians)
            if newton_raphason is None:
                E = update_exposure_gradient(E, global_gradients, t_min)
            else:
                E = update_exposure_NR(E, global_gradients, topt, tedge, newton_raphason)
        
        # Ensure non-negativity
        E = np.maximum(E, 0)
        
        # Compute loss and check convergence
        loss = -poisson.logpmf(M, (E @ S) * O)
        mean_loss = np.mean(loss)
        
        if convergence(mean_loss, loss_hat):
            print("Convergence detected")
            if conv_iter_1 == -1:
                conv_iter_1 = step
                conv_check = 0
            else:
                conv_check += 1
                
            if conv_check == 5:
                print("Algorithm converged")
                break
        else:
            conv_iter_1 = -1
            conv_check = 0
            
        loss_hat = mean_loss
        
    if conv_check < 5:
        print(f"Maximum iterations ({n_steps}) reached without convergence")
        
    return E


def running_simulation_denovo(E: np.ndarray, M: np.ndarray, S: np.ndarray, O: np.ndarray,
                            topt: float, tedge: float, lambd: float, n_steps: int) -> np.ndarray:
    """Run optimization for de novo signature discovery.
    
    Args:
        E: Initial exposure matrix
        M: Mutation count matrix
        S: Signature matrix
        O: Weight matrix
        topt: Initial optimal step size
        tedge: Initial edge step size
        lambd: L1 regularization parameter
        n_steps: Maximum number of optimization steps
    
    Returns:
        np.ndarray: Optimized exposure matrix
    """
    for step in range(n_steps):
        if np.any(E < 0):
            E = np.maximum(E, 0)
            
        # Compute gradients and step sizes
        local_gradients = compute_local_gradients(E, M, S, O)
        hessians = compute_hessians(E, M, S)
        global_gradients = compute_global_gradient(E, local_gradients, lambd)
        topt = compute_topt_denovo(E, local_gradients, global_gradients, hessians)
        tedge = compute_t_edge_denovo(E, global_gradients)
        t_min = min_topt_tedge(topt, tedge)
        
        # Update exposures
        E = update_exposure_gradient(E, global_gradients, t_min)
        
        if topt >= tedge:
            # Gradient descent update
            E = update_exposure_gradient(E, global_gradients, t_min)
        else:
            # Try Newton-Raphson update
            newton_raphason = newton_raphson1(E, global_gradients, hessians)
            if newton_raphason is None:
                E = update_exposure_gradient(E, global_gradients, t_min)
            else:
                E = update_exposure_NR(E, global_gradients, topt, tedge, newton_raphason)
        
        # Ensure non-negativity
        E = np.maximum(E, 0)
        
    return E
