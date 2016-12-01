"""
vbgmm.pyx

simple cython wrapper for variational Gaussian mixture model in C 

"""

import cython

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# declare the interface to the C code
cdef extern void c_vbgmm_fit (double* adX, int nN, int nD, int nK, int* anAssign, int debug, int bAssign)
@cython.boundscheck(False)
@cython.wraparound(False)

def fit(np.ndarray[double, ndim=2, mode="c"] xarray not None, nClusters, debug):
    """
    fit (xarray, assign, nK, debug)

    Takes a numpy array xarray as input, fits the vbgmm using nK initial clusters

    and returns cluster assignments in assign

    param: xarray -- a 2-d numpy array of np.float64
    param: assigns -- cluster assignments must have same number of rows as xarray

    """
    cdef int nN, nD, nK, nS, bAssign
        
    nN, nD = xarray.shape[0], xarray.shape[1]

    nK = nClusters

    bAssign = 0
    
    cdef np.ndarray[int, ndim=1,mode="c"] assign = np.zeros((nN), dtype=np.intc)
    
    c_vbgmm_fit (&xarray[0,0], nN, nD, nK, &assign[0], debug, bAssign)

    return assign
