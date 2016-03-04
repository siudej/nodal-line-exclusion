"""
Use sparse Cholesky decomposition to check if matrices are positive-definite.

Use CHOLMOD from SuiteSparse wrapped with scikits.sparse.
A - (12.25 + 3/256) B should be positive-definite, assuming properly scaled
stiffness matrix A and mass matrix B.

Positivity gives lower bound for the smallest eigenvalue.

Permuted and scaled A - (12.25 + 3/256) B is saved as a matrix with integer
entries for further interval arithmetic and rational computations.
"""

from scikits.sparse.cholmod import cholesky
import scipy.sparse as ss
import numpy as np


def get_matrix(filename):
    """ Import matrix from file. """
    f = open(filename)
    v = eval(f.readline())
    f.close()
    ind = np.array(v, dtype=int)[:, :2]
    print 'Values as floats: ', np.unique(np.array(v)[:, 2])
    val = np.array(v, dtype=np.int64)[:, 2]
    print 'Values as int: ', np.unique(val)
    return ss.csc_matrix((val, (ind[:, 0], ind[:, 1])))


text = open('results/cholesky.txt', 'w')
for number in range(1, 21):
    print 'Domain {}:'.format(number)
    B = get_matrix(
        './results/matrices_Z/domains_64_0_12.25_mass_{}.txt'.format(number))
    A = get_matrix(
        './results/matrices_Z/domains_64_0_12.25_stiff_{}.txt'.format(number))
    # To get eigenvalue about 12.25 we need to:
    #    * 2883  : integer entries in A, B
    #    / 128^2 : mesh of size 93 x 128 instead of 93/128 x 1
    # We scale again to get integer entries after subtracting  12.25 + 3/256,
    positive = 128 * 128 * 256 * A - (49 * 64 + 3) * 2883 * B
    #           128^2 * A - (12.25 + 3/256) * 2883 * B
    #           everything * 256 to get integers

    # cholesky factorization
    factor = cholesky(positive)
    diag = factor.D()
    L = factor.L_D()[0]
    P = factor.P()

    print >>text, 'Number of entries: ', len(positive.data)
    # the following should be positive
    print >>text, "Smallest diagonal element in LDL': ", diag.min()
    print >>text, 'Ratio largest/smallest :  ', diag.max() / diag.min()
    print >>text, 'Number of entries in L : ', len(L.data)
    print >>text, 'Permutation : ', P
    # save permuted positive
    posP = positive[:, P].tocsr()[P].tocsc()
    f = open('results/CHOLMOD_permuted/permuted{}.txt'.format(number), 'w')
    M = posP.tocoo()
    print >>f, zip(M.row, M.col, M.data)
    f.close()

    print >>text, '\n'
text.close()
