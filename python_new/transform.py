import numpy as np


def xlog2_int(n):
    return n.bit_length()-1

def truncate(s, eps):
    """
    Find the smallest k such that sum(s[:k]**2) \geq 1-eps.
    """
    mysum = 0.0
    k=-1
    while (mysum < 1-eps):
        k += 1
        mysum += s[k]**2
    return k+1


def compress(k, psi, eps=0.00001):
    """
    Compress the state so that the reduced density matrix over the first
    k qubits is the same. (But the rest is compressed.)
    """
    n=psi.size
    mat = psi.reshape(1<<k, n//(1<<k))
    U, s, Vd = np.linalg.svd(mat, full_matrices=False)
#    print(sum(s**2))
#    print(s**2)
    n_t = truncate(s, eps)
#    print(n_t)
    out = U[:,:n_t] * s[:n_t]
    
    return out.reshape(out.size, 1)


def transform_xlogx(mat):
    """
    Args:
        mat(np.array): A two-dimensional array

    Returns:
        np.array: Let UsV^† be the SVD of mat. Returns
                  Uf(s)V^†, where f(x) = -2xlogx
    """
    U, s, Vd = np.linalg.svd(mat, full_matrices=False)
    return (U * (-2.0*s * np.log(s))) @ Vd


def apply_modular_op(k, n, psi):
    """
    Given a one-dimensional vector psi over n qubits, apply the
    xlogx transform to the first k qubits.

    Args:
        k(int): Number of qubits which will undergo the xlogx transform
        n(int): Dimension
        psi(np.array): State vector

    Returns:
        np.array: Transformed state vector
    """
#    if psi.ndim>1:
#        raise TypeError("psi must be a state vector. Shape={}, Ndim={}".format(psi.shape, psi.ndim))
#    if (1<<n != psi.size):
#        raise ValueError("2**n is not equal to the dimension of psi.")

    return transform_xlogx(psi.reshape(1<<k, n//(1<<k))).reshape(n, 1)


def transform_AB(a, b, c, psi):
    """
    Apply xlogx transform to subsystem AB.

    Args:
        a(int): # of qubits in A
        b(int): # of qubits in B
        c(int): # of qubits in C
    
    Returns:
        np.array: Transformed state vector
    """
    n = psi.size
    return apply_modular_op(a+b, n, psi)


def transform_BC(a, b, c, psi):
    """
    Apply xlogx transform to subsystem BC.

    Args:
        a(int): # of qubits in A
        b(int): # of qubits in B
        c(int): # of qubits in C
    
    Returns:
        np.array: Transformed state vector
    """
    n = psi.size
    return np.transpose(apply_modular_op(b+c,n,np.transpose(psi.reshape(1<<a, n//(1<<a))).reshape(n,1)).reshape(n//(1<<a),1<<a)).reshape(n, 1)
