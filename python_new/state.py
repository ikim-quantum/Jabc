import numpy as np
import amplitude as amp


def normalize(psi):
    """
    Note: While we assume that psi is a state vector, this works
    for matrices as well if we use Hilbert-Schmidt norm.

    Args:
        psi(np.array): Unnormalized state vector

    Returns:
        np.array: Normalized state vector
    """
    norm = np.sum(np.abs(psi**2))
    return psi / np.sqrt(norm)


def rand(n):
    """
    Generate a random n-qubit state
    
    Args:
        n(int): Number of qubits

    Returns:
        np.array: Normalized random n-qubit state
    """
    psi = np.random.rand(1<<n, 1) + 1j * np.random.rand(1<<n, 1)
    return normalize(psi)
    

def amplitude(zs, x):
    """
    Given a set of coordinates, generate the wavefunction amplitude.
    (alpha=1/2)

    Args:
        zs(list(complex)): List of complex numbers

    Returns:
        complex: Amplitude
    """
    N = len(zs)
    if N%2==1:
        return 0
    elif bin(x).count("1")!=N/2:
        return 0

    znm = [zs[n] - zs[m] for n in range(N) for m in range(n+1, N) if ((x >> n & 1) == (x >> m & 1))]

    return np.prod(znm)


def amplitude_fast(zs, x):
    return amp.amplitude(zs, x)


def ansatz(zs):
    """
    Given a set of coordinates, generate the wavefunction.
    (alpha=1/2)

    Args:
        zs(list(complex)): List of complex numbers

    Returns:
        np.array: Wavefunction
    """
    N = len(zs)
    psi = np.array([amplitude_fast(zs, x) for x in range(1<<N)])

    return normalize(psi)

    
