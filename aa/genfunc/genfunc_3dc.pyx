#cython: boundscheck=False
#cython: wraparound=False

import numpy as np
cimport numpy as np
import time
from numpy.linalg import solve
from itertools import product
from scipy.linalg import cho_factor, cho_solve

def solver_cy(np.ndarray[double, ndim=2] AA, int N_max, sym=np.array([False,False,False])):
    t=time.time()
    cdef np.ndarray[double, ndim=1] NNX = np.arange(-N_max,N_max+1, 1+1*sym[0], dtype=np.double)
    cdef np.ndarray[double, ndim=1] NNY = np.arange(-N_max,N_max+1, 1+1*sym[1], dtype=np.double)
    cdef np.ndarray[double, ndim=1] NNZ = np.arange(-N_max,N_max+1, 1+1*sym[2], dtype=np.double)
    cdef double i,j,k
    cdef np.ndarray[double, ndim=2] n_vectors = np.array([[i,j,k] for (i,j,k) in product(NNX,NNY,NNZ)
                          if(not(i==0 and j==0 and k==0)    # exclude zero vector
                          and (k>0                          # northern hemisphere 
                          or (k==0 and j>0)                 # half of x-y plane
                          or (k==0 and j==0 and i>0))       # half of x axis
                          and np.sqrt(i*i+j*j+k*k)<=N_max)], dtype = np.double) # inside sphere

    cdef unsigned int nvLen = len(n_vectors)
    cdef unsigned int n = nvLen+3
    cdef np.ndarray[double] b = np.zeros(shape=(n, ), dtype=np.double)
    cdef np.ndarray[double, ndim=2] a = np.zeros(shape=(n,n), dtype=np.double)

    cdef int AALen = len(AA)
    a[:3,:3]=AALen*np.identity(3)

    cdef unsigned int n_index, aa_index, r, s, nv_index, K
    for aa_index in xrange(AALen):
        for r in xrange(3):
            for s in xrange(3,n):
                a[r,s]+=2.*n_vectors[s-3][r]*np.cos(np.dot(n_vectors[s-3],AA[aa_index][3:]))
        for r in xrange(3,n):
            for s in xrange(r,n):
                a[r,s]+=4.*np.dot(n_vectors[r-3],n_vectors[s-3])*(np.cos(np.dot(n_vectors[s-3],AA[aa_index][3:]))*np.cos(np.dot(n_vectors[r-3],AA[aa_index][3:])))

        for K in xrange(3):
            b[K]+=AA[aa_index][K]
        for n_index in xrange(nvLen):
            b[n_index+3]+=2.*np.dot(n_vectors[n_index],AA[aa_index][:3])*np.cos(np.dot(n_vectors[n_index],AA[aa_index][3:]))

    a = (a+a.T)
    for p in range(n):
        a[p,p]/=2.

    print(time.time()-t,np.shape(a))
    return np.array(cho_solve(cho_factor(a),b)), n_vectors

from itertools import izip

cdef np.ndarray[double, ndim=2] unroll_angles(np.ndarray[double, ndim=2] A):
    cdef np.ndarray[double, ndim=2] AA = np.array([0,0,0],dtype=np.double)
    cdef np.ndarray[double, ndim=2] P = np.zeros(np.shape(A),dtype=np.double)
    P[0]=A[0]
    cdef unsigned int i
    cdef np.ndarray[double] n
    for i in xrange(1,len(A)):
        n = n+(A[i]<A[i-1])*np.ones(3)*2.*np.pi
        P[i] = A[i]+n
    return P

def angle_solver_cy(np.ndarray[double, ndim=2] AA, np.ndarray[double, ndim=1] timeseries, int N_max):

    cdef np.ndarray[double, ndim=1] NN = np.arange(-N_max,N_max+1, 2,dtype=np.double)
    cdef double i,j,k
    cdef np.ndarray[double, ndim=2] n_vectors = np.array([[i,j,k] for (i,j,k) in product(NN,NN,NN)
                          if(not(i==0 and j==0 and k==0)    # exclude zero vector
                          and (k>0                          # northern hemisphere 
                          or (k==0 and j>0)                 # half of x-y plane
                          or (k==0 and j==0 and i>0))       # half of x axis
                          and np.sqrt(i*i+j*j+k*k)<=N_max)], dtype = np.double) # inside sphere

    cdef int n = 3*len(n_vectors)+6

    cdef np.ndarray[double] b = np.zeros(shape=(n), dtype=np.double)
    cdef np.ndarray[double, ndim=2] a = np.zeros(shape=(n,n), dtype=np.double)

    cdef int K, J
    for K in xrange(3):
        a[K,K]=len(AA)
        a[K,K+3]=np.sum(timeseries)
        a[K+3,K+3]=np.sum(timeseries*timeseries)

    cdef np.ndarray[double] aa,l
    cdef double TT

    for aa,TT in izip(AA,timeseries):
        for K,l in enumerate(n_vectors):
            a[6+3*K:9+3*K,:3]+=-2.*np.sin(np.dot(l,aa[3:]))*np.identity(3)
            a[6+3*K:9+3*K,3:6]+=-2.*TT*np.sin(np.dot(l,aa[3:]))*np.identity(3)
            for p,m in enumerate(n_vectors[:K+1]):
                a[6+3*K:9+3*K,6+3*p:9+3*p]+=4.*np.sin(np.dot(l,aa[3:]))*np.sin(np.dot(m,aa[3:]))*np.identity(3)
        for J in xrange(3):
            b[J]+=aa[3+J]
            b[3+J]+=TT*aa[3+J]
            for K,l in enumerate(n_vectors):
                b[6+3*K+J]+=-2.*aa[3+J]*np.sin(np.dot(l,aa[3:]))

    a = (a+a.T)
    for K in xrange(n):
        a[K,K]/=2.

    return np.array(cho_solve(cho_factor(a),b))
