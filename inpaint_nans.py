import numpy as np
import numpy.matlib as matlib
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

def inpaint_nans(A,method=0):

    if method not in [0,1,2,3,4,5]:
        raise Exception('If supplied, method must be one of:{0,1,2,3,4,5}.')

    if method != 2:
        raise Exception(f'Method {method} has not yet been implemented')

    n, m = A.shape
    A=A.flatten('F')
    nm = n * m
    k = np.isnan(A)

    nan_list = np.where(k)[0]
    known_list = np.where(~k)[0]

    nan_count = len(nan_list)

    nr = np.zeros([nan_count])
    nc = np.zeros([nan_count])

    for i in np.arange(nan_count):
        nr[i] = np.floor_divide(nan_list[i], m)
        nc[i] = nan_list[i] % m

    nan_list=np.column_stack((nan_list, nr, nc))

    if method == 2:
        # Direct solve for del^2 BVP across holes

        # generate sparse array with second partials on row
        # variable for each NaN element, only for those nodes
        # which have a row index > 1 or < n

        # is it 1-d or 2-d?
        if m == 1 or n == 1:
            # really just a 1-d case
            raise Exception('Method 2 has problems for vector input. Please use another method.')
        L = np.where(np.logical_and(nan_list[:,1] > 0, nan_list[:,1] < n-1))[0]
        nl = len(L)
        if nl > 0:
            I = matlib.repmat(nan_list[L,0].reshape(-1,1),1,3).astype(int)
            J = (matlib.repmat(nan_list[L,0].reshape(-1,1),1,3)+matlib.repmat([-1, 0, 1],nl,1)).astype(int)
            S = matlib.repmat([1, -2, 1],nl,1)
            fda = csr_matrix((S.flatten(),(I.flatten(),J.flatten())),shape=(n*m,n*m))
        else:
            fda = csr_matrix((n*m,n*m))

        # 2nd partials on column index
        L = np.where(np.logical_and(nan_list[:,2] > 0, nan_list[:,2] < m-1))[0]
        nl = len(L)
        if nl > 0:
            I = matlib.repmat(nan_list[L,0].reshape(-1,1),1,3).astype(int)
            J = matlib.repmat(nan_list[L,0].reshape(-1,1), 1,3) + matlib.repmat([-n,0,n],nl,1).astype(int)
            S = matlib.repmat([1, -2, 1],nl,1)
            fda = fda + csr_matrix((S.flatten(), (I.flatten(),J.flatten())), shape=(n*m,n*m))

        # fix boundary
        if 1 in nan_list[:,0]:
            fda[0,[0,1,n]] = np.array([-2,1,1])
        if n in nan_list[:,0]:
            fda[n-1,[n-1,n-2,(n-1)+(n-1)]] = np.array([-2,1,1])
        if nm-n+1 in nan_list[:,0]:
            fda[nm-n, [nm-n,nm-n+1,nm-n-1]] = np.array([-2,1,1])
        if nm in nan_list[:,0]:
            fda[nm, nm,nm-1,nm-n] = np.array([-2,1,1])

        # eliminate knowns
        rhs = -fda[:,known_list] * A[known_list]

        # and solve...
        B = A
        k = nan_list[:,0].astype(int)
        a = fda[k,:].tocsc()[:,k]
        b = rhs[k]
        x=spsolve(a,b)
        B[k] = x
        B = B.reshape(n,m).T
        return B