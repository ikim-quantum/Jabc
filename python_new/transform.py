import numpy as np


def xlog2_int(n):
    return n.bit_length()-1


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
        n(int): Total number of qubits
        psi(np.array): State vector

    Returns:
        np.array: Transformed state vector
    """
    if min(psi.shape)>1:
        raise TypeError("psi must be a state vector.")
    elif (1<<n != psi.size):
        raise ValueError("2**n is not equal to the dimension of psi.")

    return transform_xlogx(psi.reshape(1<<k, 1<<(n-k))).reshape(1<<n, 1)


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
    n = xlog2_int(psi.size)
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
    n = xlog2_int(psi.size)
    return np.transpose(apply_modular_op(b+c,n,np.transpose(psi.reshape(1<<a, 1<<(n-a))).reshape(1<<n,1)).reshape(1<<(n-a),1<<a)).reshape(1<<n, 1)
