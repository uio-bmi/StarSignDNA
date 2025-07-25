#!/usr/bin/env python3
"""
Main module for de novo mutational signature extraction.

This module implements the core algorithms for extracting mutational signatures
from mutation count matrices using gradient-based optimization with regularization.
The algorithm alternates between updating exposures and signatures using
gradient descent and Newton-Raphson methods.
"""

import numpy as np
from numpy.random.mtrand import poisson
from scipy.stats import poisson, entropy


def compute_local_gradients(E, M, S, O, lambd):
    """
    Compute local gradients for exposure matrix optimization.
    
    Args:
        E (np.ndarray): Current exposure matrix (n_samples x n_signatures)
        M (np.ndarray): Mutation count matrix (n_samples x n_mutations)
        S (np.ndarray): Signature matrix (n_signatures x n_mutations)
        O (np.ndarray): Opportunity matrix (n_samples x n_mutations)
        lambd (float): Regularization parameter
    
    Returns:
        tuple: (local_gradients, matrice_lambda)
            - local_gradients: Gradient matrix for exposures
            - matrice_lambda: Sample-specific regularization parameters
    """
    n_samples, n_signatures, n_mutations = (E.shape[0], S.shape[0], M.shape[1])
    local_gradients = np.empty_like(E)
    
    # Compute total mutations per sample
    matrice_sum = M.sum(axis=1, keepdims=True)
    matrice_lambda = np.empty_like(matrice_sum).astype(float)
    
    # Set regularization based on mutation count (no regularization for samples with <100 mutations)
    for idx, sum_ in np.ndenumerate(matrice_sum):
        if sum_ < 100:
            matrice_lambda[idx] = 0
        else:
            matrice_lambda[idx] = lambd
    
    # Compute gradients for each sample and signature
    for i in range(n_samples):
        for r in range(n_signatures):
            numerator = M[i] * S[r]
            denumerator_sum = np.array([E[i] @ S[:, k] for k in range(n_mutations)])
            denumerator_sum_c = denumerator_sum + 0.000001  # Add small constant for numerical stability
            local_gradients[i, r] = np.sum((numerator / denumerator_sum_c) - O[i] * S[r])
    
    return local_gradients, matrice_lambda


def compute_hessians(E, M, S):
    """
    Compute Hessian matrices for Newton-Raphson optimization.
    
    Args:
        E (np.ndarray): Current exposure matrix (n_samples x n_signatures)
        M (np.ndarray): Mutation count matrix (n_samples x n_mutations)
        S (np.ndarray): Signature matrix (n_signatures x n_mutations)
    
    Returns:
        np.ndarray: Hessian matrices for each sample
    """
    denominatior = (E @ S) ** 2 + 0.000001  # Add small constant for numerical stability
    numerator = M[:, None, None, :] * S[None, :, None, :] * S[None, None, :, :]
    res = numerator / denominatior[:, None, None, :]
    hessians = -res.sum(axis=-1)
    return hessians


def compute_global_gradient(E, local_gradients, matrice_lambda):
    """
    Compute global gradients with regularization.
    
    Args:
        E (np.ndarray): Current exposure matrix
        local_gradients (np.ndarray): Local gradients
        matrice_lambda (np.ndarray or float): Regularization parameters
    
    Returns:
        np.ndarray: Global gradients with regularization
    """
    out = np.empty_like(E)

    # Handle scalar regularization parameter
    if np.isscalar(matrice_lambda):
        for i in range(E.shape[0]):
            for j in range(E.shape[1]):
                if E[i, j] == 0 and np.all(np.abs(local_gradients[i, :]) > matrice_lambda):
                    out[i, j] = local_gradients[i, j] - matrice_lambda * np.sign(local_gradients[i, j])
                elif E[i, j] > 0:
                    out[i, j] = local_gradients[i, j] - matrice_lambda * np.sign(E[i, j])
                else:
                    out[i, j] = 0
    else:
        # Handle sample-specific regularization parameters
        for i in range(E.shape[0]):
            for j in range(E.shape[1]):
                if E[i, j] == 0 and np.all(np.abs(local_gradients[i, :]) > matrice_lambda[i]):
                    out[i, j] = local_gradients[i, j] - matrice_lambda[i] * np.sign(local_gradients[i, j])
                elif E[i, j] > 0:
                    out[i, j] = local_gradients[i, j] - matrice_lambda[i] * np.sign(E[i, j])
                else:
                    out[i, j] = 0

    return out


def compute_topt(E, local_gradients, global_gradients, hessians):
    """
    Compute optimal step size for gradient descent.
    
    Args:
        E (np.ndarray): Current exposure matrix
        local_gradients (np.ndarray): Local gradients
        global_gradients (np.ndarray): Global gradients
        hessians (np.ndarray): Hessian matrices
    
    Returns:
        float: Optimal step size
    """
    numerator = np.linalg.norm(global_gradients, ord=None, axis=None, keepdims=False)
    gg_vectors = (gg[:, None] for gg in global_gradients)
    denominatior = sum([gg.T @ hessians @ gg for gg, hessians in zip(gg_vectors, hessians)])
    topt = -(numerator / denominatior)
    return topt


def compute_t_edge(E, global_gradients):
    """
    Compute maximum step size to maintain non-negativity.
    
    Args:
        E (np.ndarray): Current exposure matrix
        global_gradients (np.ndarray): Global gradients
    
    Returns:
        float: Maximum step size
    """
    global_gradients = global_gradients.flatten()
    E = E.flatten()
    
    # Remove elements where gradient is zero
    ind = np.where(global_gradients == 0)
    E_Conv = np.delete(E, ind[0])
    global_gradients_conv = np.delete(global_gradients, ind[0])
    
    # Find where signs are opposite (potential boundary crossing)
    mask = np.sign(E_Conv) == -np.sign(global_gradients_conv)
    mask &= (np.sign(E_Conv) != 0)
    
    if not np.any(mask):
        return np.inf
    
    assert np.all(global_gradients_conv != 0)
    return np.min(-(E_Conv / global_gradients_conv)[mask])


def min_topt_tedge(topt, tedge):
    """
    Select the minimum of optimal and edge step sizes.
    
    Args:
        topt (float): Optimal step size
        tedge (float): Edge step size
    
    Returns:
        float: Minimum step size
    """
    topt_tedge = np.minimum(float(topt), float(tedge))
    return topt_tedge


def update_exposure_gradient(E, global_gradients, topt_tedge):
    """
    Update exposures using gradient descent.
    
    Args:
        E (np.ndarray): Current exposure matrix
        global_gradients (np.ndarray): Global gradients
        topt_tedge (float): Step size
    
    Returns:
        np.ndarray: Updated exposure matrix
    """
    if np.any(E < 0):
        E = np.maximum(E, 0)
    return E + topt_tedge * global_gradients


def newton_raphson1(E, global_gradients, hessians):
    """
    Apply Newton-Raphson optimization step.
    
    Args:
        E (np.ndarray): Current exposure matrix
        global_gradients (np.ndarray): Global gradients
        hessians (np.ndarray): Hessian matrices
    
    Returns:
        np.ndarray or None: Updated exposure matrix or None if singular
    """
    nr = []
    H = hessians
    active_mask = (E != 0)
    
    assert np.all(E >= 0), "Error: E matrix element is not greater than zero."
    
    active_hessians = []
    for E_row, gradient_row, hessian in zip(E, global_gradients, H):
        non_active = ((E_row == 0) | (gradient_row == 0))
        active_mask = ~non_active
        active_gradients = gradient_row[active_mask]
        active_hessian = hessian[active_mask][:, active_mask]
        new_row = E_row.copy()
        det = np.linalg.det(active_hessian)
        active_hessians.append(det ** -1)
        
        if det < 10e-10:
            return None
            
        new_row[active_mask] = E_row[active_mask] - np.linalg.inv(active_hessian) @ active_gradients
        nr.append(new_row)
    
    v1 = np.array(nr)
    return v1


def update_exposure_NR(E, global_gradients, topt, tedge, new_E):
    """
    Update exposures using Newton-Raphson method.
    
    Args:
        E (np.ndarray): Current exposure matrix
        global_gradients (np.ndarray): Global gradients
        topt (float): Optimal step size
        tedge (float): Edge step size
        new_E (np.ndarray): Newton-Raphson update
    
    Returns:
        np.ndarray: Updated exposure matrix
    """
    if np.any(E < 0):
        E = np.maximum(E, 0)
    return np.where(np.sign(new_E) == np.sign(E),
                    new_E,
                    E + topt * global_gradients)


def convergence(E, E_hat, tol=10e-6):
    """
    Check convergence between two exposure matrices.
    
    Args:
        E (np.ndarray): Current exposure matrix
        E_hat (np.ndarray): Previous exposure matrix
        tol (float): Convergence tolerance
    
    Returns:
        bool: True if converged, False otherwise
    """
    conv = np.abs((E_hat - E) / E)
    if conv < tol:
        return True
    else:
        return False


def mean_exposure(E):
    """
    Compute mean exposure across all samples and signatures.
    
    Args:
        E (np.ndarray): Exposure matrix
    
    Returns:
        float: Mean exposure
    """
    m = np.mean(E)
    return m


def Frobinous(M, S, E, O):
    """
    Compute Frobenius norm of reconstruction error.
    
    Args:
        M (np.ndarray): Original mutation matrix
        S (np.ndarray): Signature matrix
        E (np.ndarray): Exposure matrix
        O (np.ndarray): Opportunity matrix
    
    Returns:
        float: Frobenius norm of reconstruction error
    """
    from numpy import linalg as LA
    fibo1 = (E @ S) * O
    fibo = LA.norm(M - fibo1, ord=2)
    return fibo


def Frobinous_reconstuct(M, S, E, O):
    """
    Compute Frobenius norm of reconstruction error with count scaling.
    
    Args:
        M (np.ndarray): Original mutation matrix
        S (np.ndarray): Signature matrix
        E (np.ndarray): Exposure matrix
        O (np.ndarray): Opportunity matrix
    
    Returns:
        tuple: (frobenius_norm, reconstructed_matrix)
    """
    from numpy import linalg as LA
    fibo1 = (E @ S) * O
    M_count = M.sum(axis=1)
    fibo1 *= M_count.reshape(-1, 1)
    M_hat = fibo1
    fibo = LA.norm(M - M_hat, ord=2)
    return fibo, M_hat


def sparsity(E):
    """
    Check if exposure matrix is sparse (>=80% zeros).
    
    Args:
        E (np.ndarray): Exposure matrix
    
    Returns:
        bool: True if sparse, False otherwise
    """
    sparsity = 1.0 - (np.count_nonzero(E) / float(E.size))
    if sparsity >= 0.8:
        return True
    else:
        return False


def mse(E, E_hat):
    """
    Compute mean squared error between two matrices.
    
    Args:
        E (np.ndarray): First matrix
        E_hat (np.ndarray): Second matrix
    
    Returns:
        float: Mean squared error
    """
    mse_error = np.square(np.subtract(E, E_hat)).mean()
    return mse_error


def cos_sim_matrix(matrix1, matrix2):
    """
    Compute cosine similarity matrix and optimal assignment.
    
    Args:
        matrix1 (np.ndarray): First signature matrix
        matrix2 (np.ndarray): Second signature matrix
    
    Returns:
        tuple: (similarity_matrix, assignment_string, matched_matrix)
    """
    import scipy.spatial as sp
    import pandas as pd
    import string
    from scipy.optimize import linear_sum_assignment
    
    # Compute cosine similarity matrix
    cosine = 1 - sp.distance.cdist(matrix1, matrix2, 'cosine')
    m = 1 - sp.distance.cdist(matrix1, matrix2, 'cosine')
    
    # Find optimal assignment
    match_coresp = linear_sum_assignment(m, maximize=True)
    print(match_coresp)
    
    match_coresp_1 = "Optimal assignment: \n"
    select_row = []
    for i, j in zip(match_coresp[0], match_coresp[1]):
        match_coresp_1 = match_coresp_1 + str(i) + " => " + str(j) + ", "
        select_row.append(j)
    
    # Extract matched signatures
    cosine_transpose = np.transpose(cosine)
    cosine_extract = cosine_transpose[select_row]
    
    # Create labels for de novo signatures
    rows_names = []
    index_names = matrix2.index.values.tolist()
    for index in select_row:
        rows_names.append(index_names[index])
    
    alphabet = list(string.ascii_uppercase)
    cosine_extract_columns = ['Denovo ' + alphabet[k] for k in range(cosine_extract.shape[1])]
    
    # Create output dataframes
    df = pd.DataFrame(cosine_extract, columns=cosine_extract_columns)
    df[" "] = rows_names
    df = df.set_index([" "])
    
    cosine_index = pd.DataFrame(cosine_transpose, columns=cosine_extract_columns)
    cosine_index[" "] = index_names
    cosine_index_detect = cosine_index.set_index([" "])
    
    return cosine_index_detect, match_coresp_1[:-2], df


def running_simulation_denovo(E, M, S, O, topt, tedge, lambd):
    """
    Main optimization loop for de novo signature extraction.
    
    Args:
        E (np.ndarray): Initial exposure matrix
        M (np.ndarray): Mutation count matrix
        S (np.ndarray): Signature matrix
        O (np.ndarray): Opportunity matrix
        topt (float): Initial optimal step size
        tedge (float): Initial edge step size
        lambd (float): Regularization parameter
    
    Returns:
        np.ndarray: Optimized exposure matrix
    """
    old_loss = np.inf
    pmf_s = []
    mse_e = 0
    loss = 0
    conv_check = 0
    conv_iter_1 = 0
    mse_hat = mse_e
    
    # Main optimization loop
    for step in range(25):
        mse_hat = mse_e
        loss_hat = loss
        
        # Ensure non-negativity
        if np.any(E < 0):
            E = np.maximum(E, 0)
        
        # Compute gradients and Hessians
        local_gradients, matrice_lambda = compute_local_gradients(E, M, S, O, lambd)
        hessians = compute_hessians(E, M, S)
        global_gradients = compute_global_gradient(E, local_gradients, matrice_lambda)
        
        # Compute step sizes
        topt = compute_topt(E, local_gradients, global_gradients, hessians)
        tedge = compute_t_edge(E, global_gradients)
        minimun_topt_tedge = min_topt_tedge(topt, tedge)
        
        # Update exposures
        E = update_exposure_gradient(E, global_gradients, minimun_topt_tedge)
        
        # Choose optimization method based on step size relationship
        if topt >= tedge:
            # Use gradient descent
            if np.any(E < 0):
                E = np.maximum(E, 0)
            E = update_exposure_gradient(E, global_gradients, minimun_topt_tedge)
            mse_e = Frobinous(M, S, E, O)
            loss = -poisson.logpmf(M, (E @ S) * O)
        else:
            # Try Newton-Raphson, fall back to gradient descent if singular
            if np.any(E < 0):
                E = np.maximum(E, 0)
            newton_raphason = newton_raphson1(E, global_gradients, hessians)
            
            if newton_raphason is None:
                # Fall back to gradient descent
                if np.any(E < 0):
                    E = np.maximum(E, 0)
                E = update_exposure_gradient(E, global_gradients, minimun_topt_tedge)
            else:
                # Use Newton-Raphson update
                if np.any(E < 0):
                    E = np.maximum(E, 0)
                E = update_exposure_NR(E, global_gradients, topt, tedge, newton_raphason)
                mse_e = Frobinous(M, S, E, O)
                loss = -poisson.logpmf(M, (E @ S) * O)
        
        # Ensure non-negativity at the end
        if np.any(E < 0):
            E = np.maximum(E, 0)
    
    return E
