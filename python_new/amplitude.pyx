import numpy as np

cdef extern int __builtin_popcount(unsigned int)

def amplitude(complex[:] zs, int x):
    cdef int N = zs.size        
    cdef int n, m
    cdef complex total = 1.0

    for n in range(N):
        for m in range(n+1, N):
            if ((x>> n & 1) == (x >> m & 1)):
                total *= (zs[n] - zs[m])

    return total


def ansatz(complex[:] zs):
    cdef int N = zs.size
    out = np.zeros(1<<N, dtype=np.cdouble)

    if N%2==1:
        return out
    else:
        for x in range(1<<N):
            if __builtin_popcount(x)==N/2:
                out[x] = amplitude(zs, x)
        return out