from state import ansatz
from transform import transform_AB, transform_BC, compress
import numpy as np
import time


def wavefunction(zs_A, zs_B, zs_C, zs_D):
    zs = np.array(list(zs_A) + list(zs_B) + list(zs_C) + list(zs_D))
    
    print("Generating ansatz for {} qubits.".format(len(zs)))
    return ansatz(zs)


def Jabc(zs_A, zs_B, zs_C, zs_D, verbose=True):
    if verbose:
        start = time.time()
    psi = wavefunction(zs_A, zs_B, zs_C, zs_D)
    if verbose:
        end = time.time()
        print("Wavefunction generated. Time elapsed={} seconds".format(end-start))
    a, b, c, d= len(zs_A), len(zs_B), len(zs_C), len(zs_D)

    if verbose:
        print("AB transform")
        start=time.time()
    psi_AB = transform_AB(a, b, c, psi)
    if verbose:
        end = time.time()
        print("Done. Time elapsed={} seconds".format(end-start))

    if verbose:
        print("BC transform")
        start=time.time()
    psi = transform_BC(a, b, c, psi)
    if verbose:
        end = time.time()
        print("Done. Time elapsed={} seconds".format(end-start))

    prod = psi_AB.T.conj() @ psi
    return -2.0 *prod[0][0].imag


def Jabc2(zs_A, zs_B, zs_C, zs_D, verbose=True):
    k = len(zs_A) + len(zs_B) + len(zs_C)
    kk = k + len(zs_D)
    if verbose:
        start = time.time()
    
    psi = compress(k, wavefunction(zs_A, zs_B, zs_C, zs_D))
    if verbose:
        end = time.time()
        print("Wavefunction generated. Time elapsed={} seconds".format(end-start))

    print("Original wavefunction dimension: {} ".format(1<<kk))
    print("Compressed wavefunction dimension: {}".format(psi.size))

    a, b, c = len(zs_A), len(zs_B), len(zs_C)

    if verbose:
        print("AB transform")
        start=time.time()
    psi_AB = transform_AB(a, b, c, psi)
    if verbose:
        end = time.time()
        print("Done. Time elapsed={} seconds".format(end-start))

    if verbose:
        print("BC transform")
        start=time.time()
    psi = transform_BC(a, b, c, psi)
    if verbose:
        end = time.time()
        print("Done. Time elapsed={} seconds".format(end-start))

    prod = psi_AB.T.conj() @ psi
    return -2.0 *prod[0][0].imag
    
