import numpy as np


# function to compute the local gradient
from numpy.random.mtrand import poisson
from scipy.stats import poisson, entropy

def compute_local_gradients(E, M, S, O,lambd):
    n_samples, n_signatures, n_mutations = (E.shape[0], S.shape[0], M.shape[1])
    local_gradients = np.empty_like(E)
    matrice_sum = M.sum(axis=1, keepdims=True)
    matrice_lambda = np.empty_like(matrice_sum).astype(float)
    for idx, sum_ in np.ndenumerate(matrice_sum):
        if sum_ < 100:
            matrice_lambda[idx] = 0
        else:
            matrice_lambda[idx] = lambd
    for i in range(n_samples):
        for r in range(n_signatures):
            numerator = M[i] * S[r]
            denumerator_sum = np.array([E[i] @ S[:, k] for k in range(n_mutations)])
            denumerator_sum_c = denumerator_sum    + 0.000001
            local_gradients[i, r] = np.sum((numerator / denumerator_sum_c) - O[i] * S[r])
    return local_gradients, matrice_lambda

# function to calculate the hessian matrix
def compute_hessians(E, M, S):
    denominatior = (E @ S) ** 2 + 0.000001
    numerator = M[:, None, None, :] * S[None, :, None, :] * S[None, None, :, :]
    res = numerator / denominatior[:, None, None, :]
    # print(res)
    hessians = -res.sum(axis=-1)
    # return -res.sum(axis=-1)
    return hessians


# function to compute the global gradient
def compute_global_gradient(E, local_gradients, matrice_lambda):
    out = np.empty_like(E)

    # Check if matrice_lambda is a scalar
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
        for i in range(E.shape[0]):
            for j in range(E.shape[1]):
                if E[i, j] == 0 and np.all(np.abs(local_gradients[i, :]) > matrice_lambda[i]):
                    out[i, j] = local_gradients[i, j] - matrice_lambda[i] * np.sign(local_gradients[i, j])
                elif E[i, j] > 0:
                    out[i, j] = local_gradients[i, j] - matrice_lambda[i] * np.sign(E[i, j])
                else:
                    out[i, j] = 0

    return out

# function to compute the step-size
def compute_topt(E, local_gradients, global_gradients, hessians):
    # print("local",hessians)
    numerator = np.linalg.norm(global_gradients, ord=None, axis=None, keepdims=False)
    # numerator = np.sum(global_gradients * local_gradients)
    # print("numerator", numerator)
    gg_vectors = (gg[:, None] for gg in global_gradients)
    denominatior = sum([gg.T @ hessians @ gg for gg, hessians in zip(gg_vectors, hessians)])
    topt = - (numerator / denominatior)
    return topt


# function to compute the maximimum step-size value i.e the maximum step-size

def compute_t_edge(E, global_gradients):
    global_gradients = global_gradients.flatten()
    E = E.flatten()
    ind = np.where(global_gradients == 0)
    E_Conv = np.delete(E, ind[0])
    global_gradients_conv = np.delete(global_gradients, ind[0])
    mask = np.sign(E_Conv) == -np.sign(global_gradients_conv)
    mask &= (np.sign(E_Conv) != 0)
    if not np.any(mask):
        return np.inf
    assert np.all(global_gradients_conv != 0)
    return np.min(-(E_Conv / global_gradients_conv)[mask])


# function to select the minimum step size
def min_topt_tedge(topt, tedge):
    topt_tedge = np.minimum(float(topt), float(tedge))
    return topt_tedge


# function to update gradient Ascent algorithm
def update_exposure_gradient(E, global_gradients, topt_tedge):
    if (np.any(E < 0)):
        E = np.maximum(E, 0)
    return E + topt_tedge * global_gradients


# function to compute Newton_Raphson
def newton_raphson1(E, global_gradients, hessians):
    nr = []
    v1 = []
    H = hessians
    active_mask = (E != 0)
    # print(E)
    assert np.all(E >= 0), "Error: E matrix element is not  greater than zero."
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


# function to update newton optimization
def update_exposure_NR(E, global_gradients, topt, tedge, new_E):
    if (np.any(E < 0)):
        E = np.maximum(E, 0)
    # assert topt <= tedge
    return np.where(np.sign(new_E) == np.sign(E),
                    new_E,
                    E + topt * global_gradients)


def convergence(E, E_hat, tol=10e-6):
    conv = []
    conv = np.abs((E_hat - E) / E)
    if conv < tol:
        return True
    else:
        return False

def mean_exposure(E):
    m = []
    m = np.mean(E)
    return m


# function to compute the frobinous norm
def Frobinous(M, S, E, O):
    from numpy import linalg as LA
    fibo = []
    fibo1 = []
    fibo1 = (E @ S) * O
    fibo = LA.norm(M - fibo1, ord=2)
    return fibo


# function to compute the frobinous norm
def Frobinous_reconstuct(M, S, E, O):
    from numpy import linalg as LA
    fibo = []
    fibo1 = []
    M_count = []
    M_hat = []
    fibo1 = (E @ S) * O
    M_count = M.sum(axis=1)
    fibo1 *= M_count.reshape(-1, 1)
    M_hat = fibo1
    fibo = LA.norm(M - M_hat, ord=2)
    return fibo, M_hat

#calculation of sparsity
def sparsity(E):
    sparsity = 1.0 - (np.count_nonzero(E) / float(E.size))
    # print("The sparsity is :", sparsity)
    if sparsity >= 0.8:
        return True
    else:
        return False

def mse(E, E_hat):
    mse_error = []
    # from sklearn.metrics import mean_squared_error
    mse_error = np.square(np.subtract(E, E_hat)).mean()
    return mse_error


# function to compute the cosine similarity function
def cos_sim_matrix(matrix1, matrix2):
    import scipy.spatial as sp
    import pandas as pd
    import string
    from scipy.optimize import linear_sum_assignment
    index_names = []
    cosine = 1 - sp.distance.cdist(matrix1, matrix2, 'cosine')
    m = 1 - sp.distance.cdist(matrix1, matrix2, 'cosine')
    match_coresp = linear_sum_assignment(m, maximize=True)
    print(match_coresp)
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

#
def running_simulation_denovo(E, M, S, O, topt, tedge, lambd):
    old_loss = np.inf
    pmf_s = []
    mse_e = 0
    loss = 0
    conv_check = 0
    conv_iter_1 = 0
    mse_hat = mse_e
    for step in range(25):
        mse_hat = mse_e
        loss_hat = loss
        # print("Gradient step is:", step)
        # E_hat = E
        if np.any(E < 0):
            E = np.maximum(E, 0)
        local_gradients, matrice_lambda = compute_local_gradients(E, M, S, O, lambd)
        hessians = compute_hessians(E, M, S)
        global_gradients = compute_global_gradient(E, local_gradients, matrice_lambda)
        topt = compute_topt(E, local_gradients, global_gradients, hessians)
        tedge = compute_t_edge(E, global_gradients)
        minimun_topt_tedge = min_topt_tedge(topt, tedge)
        E = update_exposure_gradient(E, global_gradients, minimun_topt_tedge)
        if topt >= tedge:
            if np.any(E < 0):
                E = np.maximum(E, 0)
            E = update_exposure_gradient(E, global_gradients, minimun_topt_tedge)
            mse_e = Frobinous(M, S, E, O)
            loss = -poisson.logpmf(M, (E @ S) * O)
        else:
            if np.any(E < 0):
                E = np.maximum(E, 0)
            newton_raphason = newton_raphson1(E, global_gradients, hessians)
            if newton_raphason is None:
                if np.any(E < 0):
                    E = np.maximum(E, 0)
                E = update_exposure_gradient(E, global_gradients, minimun_topt_tedge)
            else:
                if np.any(E < 0):
                    E = np.maximum(E, 0)
                E = update_exposure_NR(E, global_gradients, topt, tedge, newton_raphason)
                mse_e = Frobinous(M, S, E, O)
                loss = -poisson.logpmf(M, (E @ S) * O)
        if np.any(E < 0):
            E = np.maximum(E, 0)
    return E
