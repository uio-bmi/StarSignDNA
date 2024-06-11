import numpy as np


# function to compute the local gradient
from numpy.random.mtrand import poisson
from scipy.stats import poisson, entropy

def compute_local_gradients(E, M, S, O):
    n_samples, n_signatures, n_mutations = (E.shape[0], S.shape[0], M.shape[1])
    local_gradients = np.empty_like(E)
    for i in range(n_samples):
        for r in range(n_signatures):
            # print(type(M[i]))
            # print(type(S[r]))
            # print(M.dtype, S.dtype)
            numerator = M[i]*S[r]
            denumerator_sum = np.array([E[i] @ S[:, k] for k in range(n_mutations)])
            denumerator_sum_c = denumerator_sum +  0.000001
            local_gradients[i, r] = np.sum((numerator/denumerator_sum_c) - O[i]*S[r])

    return local_gradients


def compute_local_gradients_denovo(E, M, S, O):
    n_samples, n_signatures, n_mutations = (E.shape[0], S.shape[0], M.shape[1])
    local_gradients = np.empty_like(E)
    for i in range(n_samples):
        for r in range(n_signatures):
            numerator = M[i] * S[r]
            denumerator_sum = np.array([E[i] @ S[:, k] for k in range(n_mutations)])
            denumerator_sum_c = denumerator_sum #+ 10e-14 #np.finfo(float).eps # 0.000001
            local_gradients[i, r] = np.sum((numerator / denumerator_sum_c) - O[i] * S[r])
    return local_gradients


# function to calculate the hessian matrix
def compute_hessians(E, M, S):
    denominatior = (E @ S) ** 2 + np.finfo(float).eps # 0.000001
    numerator = M[:, None, None, :] * S[None, :, None, :] * S[None, None, :, :]
    res = numerator / denominatior[:, None, None, :]
    # print(res)
    #print("hessian", denominatior)
    hessians = -res.sum(axis=-1)

    # return -res.sum(axis=-1)
    return hessians



def compute_global_gradient(E, local_gradients, lambd):
    cond_a = local_gradients-lambd*np.sign(E) #change 7feb
    cond_b = local_gradients-lambd*np.sign(local_gradients) #change 7feb
    cond_c = 0
#    # print(np.where(E!=0, cond_a, np.where(np.abs(local_gradients)>lambd, cond_b, cond_c)))
    return np.where(E != 0, cond_a, np.where(np.abs(local_gradients) > lambd, cond_b, cond_c))




def compute_global_gradient_denovo(E, local_gradients, matrice_lambda):
    cond_a = local_gradients - matrice_lambda * np.sign(E)
    cond_b = local_gradients - matrice_lambda * np.sign(local_gradients)
    cond_c = 0
    return np.where(E != 0, cond_a, np.where(np.all(np.abs(local_gradients) > matrice_lambda), cond_b, cond_c))


# function to compute the step-size
def compute_topt(E, local_gradients, global_gradients, hessians):
    if (np.any(E < 0)):
        E = np.maximum(E, 0)
    # print("local",hessians)
    numerator = np.linalg.norm(global_gradients, ord= 'fro', axis=None, keepdims=False)
    ##numerator = np.sum(global_gradients * local_gradients)
    # print("numerator", numerator)
    gg_vectors = (gg[:, None] for gg in global_gradients)
    denominatior = sum([gg.T @ hessians @ gg for gg, hessians in zip(gg_vectors, hessians)])
    topt = - (numerator / denominatior)  + 10e-5
    return topt


# function to compute the maximimum step-size value i.e the maximum step-size



def compute_t_edge(E, global_gradients):
    global_gradients = global_gradients.flatten()
    E = E.flatten()
    mask = np.sign(E) == -np.sign(global_gradients)
    mask &= (np.sign(E) != 0)
    if not np.any(mask):
        return np.inf
    return np.min(-(E/global_gradients)[mask])  + 10e-2

def compute_topt_denovo(E, local_gradients, global_gradients, hessians):
    # print("local",hessians)
    #numerator = np.linalg.norm(global_gradients, ord='fro', axis=None, keepdims=False)
    numerator = np.sum(global_gradients * local_gradients)
    # print("numerator", numerator)
    gg_vectors = (gg[:, None] for gg in global_gradients)
    denominatior = sum([gg.T @ hessians @ gg for gg, hessians in zip(gg_vectors, hessians)])
    topt = - (numerator / denominatior)
    return topt


# function to compute the maximimum step-size value i.e the maximum step-size

def compute_t_edge_denovo(E, global_gradients):
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
    ##topt_tedge = np.min(topt, tedge)
    return topt_tedge


# function to update gradient Ascent algorithm
def update_exposure_gradient(E, global_gradients, topt_tedge):
    # print("topt_tedge",topt_tedge)
    if (np.any(E < 0)): ##25feb change
       E = np.maximum(E, 0) ##feb change
    return E + topt_tedge * global_gradients


# function to compute Newton_Raphson
def newton_raphson1(E, global_gradients, hessians):
    nr = []
    v1 = []
    H = hessians
    active_mask = (E != 0)
    # print(E)
    assert np.all(E >= 0), "Error: E matrix element is not  greater than zero." #comment 25jan
    active_hessians = []
    for E_row, gradient_row, hessian in zip(E, global_gradients, H):
        non_active = ((E_row == 0) | (gradient_row == 0))
        active_mask = ~non_active
        active_gradients = gradient_row[active_mask]
        active_hessian = hessian[active_mask][:, active_mask]
        #print("active hessian", active_hessian)
        #print("active hessian_shape", active_hessian.shape)
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
       E = np.maximum(E, 0) ##25jan change
       assert topt <= tedge
   return np.where(np.sign(new_E) == np.sign(E),
                    new_E,
                    E + topt * global_gradients)


def convergence(E, E_hat, tol=1e-8):
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
    fibo = LA.norm(M - fibo1, ord=None)
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
    fibo = LA.norm(M - M_hat, ord=None)
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
  #  return  cosine_index_detect

# function to run the optimisation algorithm
def running_simulation_refit(E, M, S, O, topt, tedge, lambd, n_steps):
    old_loss = np.inf
    pmf_s = []
    mse_e = 0
    min_exo_f = []
    min_expo_in = 0
    loss = 0
    # mse_hat = 0
    expose_mat = []
    mse_s = []
    expo_step = []
    expose_mat1 = []
    conv_iter_1 = 0
    conv_check = 0
    mse_old = np.inf
    loss_hat = np.inf
    minimun_topt_tedge = 0
    for step in range(n_steps):
        print("Gradient Step is:", step)
     ###   loss = -poisson.logpmf(M, (E @ S) * O)
        mse_hat = mse_e
        # loss_hat = loss
        E_hat = E
        if np.any(E < 0):  #comment 25jan
            E = np.maximum(E, 0) #comment 25jan
        ###local_gradients, matrice_lambda = compute_local_gradients(E, M, S, O, lambd)
        local_gradients = compute_local_gradients(E, M, S, O)
        #print(matrice_lambda)
        #print(local_gradients)
        #print(E)
        hessians = compute_hessians(E, M, S)
        ###global_gradients = compute_global_gradient(E, local_gradients, matrice_lambda) #threshold lambda
        global_gradients = compute_global_gradient(E, local_gradients,lambd)
        tedge = compute_t_edge(E, global_gradients)
        topt = compute_topt(E, local_gradients, global_gradients, hessians)
        #print("topt", topt)
        #print("tedge",tedge)
        minimun_topt_tedge = min_topt_tedge(topt, tedge)
        #print("Minim", minimun_topt_tedge)
        #    print("minimun_topt_tedge",minimun_topt_tedge)
        E = update_exposure_gradient(E, global_gradients, minimun_topt_tedge)
        if topt >= tedge:

            ##add24april
            local_gradients = compute_local_gradients(E, M, S, O)
            hessians = compute_hessians(E, M, S)
            global_gradients = compute_global_gradient(E, local_gradients,lambd)
            ####add24april
            #topt = tedge
            #print("HERRRRE")
            minimun_topt_tedge = min_topt_tedge(topt, tedge)
        #    print("minimun_topt_tedge",minimun_topt_tedge)
            E = update_exposure_gradient(E, global_gradients, minimun_topt_tedge)
            mse_e = Frobinous(M, S, E, O)
            loss = -poisson.logpmf(M, (E @ S) * O)
            # min_expo_in = np.min(E)
            if np.any(E < 0): #comment 25jan
                E = np.maximum(E, 0) #comment 25jan
        else:
            if np.any(E < 0): #comment 25jan
                E = np.maximum(E, 0) #comment 25jan
            newton_raphason = newton_raphson1(E, global_gradients, hessians)
            if newton_raphason is None:
                #print("NEWTON1")
                if np.any(E < 0): #comment 25jan
                    E = np.maximum(E, 0) #comment 25jan
                E = update_exposure_gradient(E, global_gradients, minimun_topt_tedge)
                mse_e = Frobinous(M, S, E, O)
                loss = -poisson.logpmf(M, (E @ S) * O)
                if np.any(E < 0): #comment 25jan
                    E = np.maximum(E, 0) #comment 25jan
            else:
                #print("NEWTON2")
                if np.any(E < 0): #comment 25jan
                    E = np.maximum(E, 0) #comment 25jan
                topt = compute_topt(E, local_gradients, global_gradients, hessians)
                tedge = compute_t_edge(E, global_gradients)
                E = update_exposure_NR(E, global_gradients, topt, tedge, newton_raphason)
                mse_e = Frobinous(M, S, E, O)
                loss = -poisson.logpmf(M, (E @ S) * O)
        if np.any(E < 0): #comment 25jan
            E = np.maximum(E, 0) #comment 25jan
        conv = convergence(np.mean(loss_hat), np.mean(loss))
        loss_hat = np.mean(loss)
        if conv == True:
            print(f"StarSign converge: {conv}")
            if conv_iter_1 == -1:
                conv_iter_1 = step
                conv_check = 0
            else:
                conv_check = conv_check + 1
        else:
            print(f" StarSign converged: {conv}")
            conv_iter_1 = -1
            conv_check = 0
        if conv_check == 5:
            print("Thanks: StarSign Algorithm converged")
            break
        mse_old = mse_e
        loss_hat = np.mean(loss)
        loss = -poisson.logpmf(M, (E @ S) * O)
        mse_e = Frobinous(M, S, E, O)
        loss[loss == np.float64("Inf")] = 0
        pmf_s.append(np.mean(loss))
        mse_s.append(np.mean(mse_e))
        expo_step.append(np.mean(E, axis=0))
    return E

def running_simulation_denovo(E, M, S, O, topt, tedge, lambd, n_steps):
    old_loss = np.inf
    pmf_s = []
    mse_e = 0
    loss = 0
    conv_check = 0
    conv_iter_1 = 0
    mse_hat = mse_e
    for step in range(n_steps):
        mse_hat = mse_e
        loss_hat = loss
        # print("Gradient step is:", step)
        # E_hat = E
        if np.any(E < 0):
            E = np.maximum(E, 0)
        local_gradients = compute_local_gradients(E, M, S, O)
      #  local_gradients, matrice_lambda = compute_local_gradients(E, M, S, O, lambd)
        hessians = compute_hessians(E, M, S)
        global_gradients = compute_global_gradient(E, local_gradients, lambd)
     ##   global_gradients = compute_global_gradient(E, local_gradients, matrice_lambda)
        topt = compute_topt_denovo(E, local_gradients, global_gradients, hessians)
        tedge = compute_t_edge_denovo(E, global_gradients)
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
