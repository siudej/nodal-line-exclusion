"""
Find interval arithmetic based LDL' decomposition for matrices.

Using sympy SparseMatrix as a container, coerced to accept mpmath mpi
intervals as entries.

LDL decomposition uses sympy generated sparsity pattern, and a modified
_LDL_sparse method from SparseMatrix.
"""

from sympy.matrices import SparseMatrix
import scipy.sparse as ss
import numpy as np
from mpmath import mpi


def LDL(mat):
    """
    Algorithm for numeric LDL factization, exploiting sparse structure.

    This function is a modification of scipy.sparse.SparseMatrix._LDL_sparse,
    allowing mpmath.mpi interval arithmetic objects as entries.

    L, D are SparseMatrix objects. However we assign values through _smat member
    to avoid type conversions to Rational.
    """
    Lrowstruc = mat.row_structure_symbolic_cholesky()
    print 'Number of entries in L: ', np.sum(map(len, Lrowstruc))
    L = SparseMatrix(mat.rows, mat.rows,
                     dict([((i, i), mpi(0)) for i in range(mat.rows)]))
    D = SparseMatrix(mat.rows, mat.cols, {})
    for i in range(len(Lrowstruc)):
        for j in Lrowstruc[i]:
            if i != j:
                L._smat[(i, j)] = mat._smat.get((i, j), mpi(0))
                summ = 0
                for p1 in Lrowstruc[i]:
                    if p1 < j:
                        for p2 in Lrowstruc[j]:
                            if p2 < j:
                                if p1 == p2:
                                    summ += L[i, p1]*L[j, p1]*D[p1, p1]
                            else:
                                break
                    else:
                        break
                L._smat[(i, j)] = L[i, j] - summ
                L._smat[(i, j)] = L[i, j] / D[j, j]

            elif i == j:
                D._smat[(i, i)] = mat._smat.get((i, i), mpi(0))
                summ = 0
                for k in Lrowstruc[i]:
                    if k < i:
                        summ += L[i, k]**2*D[k, k]
                    else:
                        break
                D._smat[(i, i)] -= summ

    return L, D


def get_matrix(filename):
    """
    Import matrix from file.

    Entries assumed to be integers.
    """
    f = open(filename)
    v = np.array(eval(f.readline()), dtype=int)
    f.close()
    ind = np.array(v, dtype=int)[:, :2]
    val = np.array(v)[:, 2]
    return ss.csc_matrix((val, (ind[:, 0], ind[:, 1])))


text = open('results/ldl_mpi.txt', 'w')
for number in range(1, 21):
    positive = get_matrix('./results/CHOLMOD_permuted/permuted{}.txt'.format(number))
    size = max(positive.indices)+1
    print >>text, 'Domain {}:'.format(number)
    print >>text, 'Number of entries: ', len(positive.data)
    # change integers into mpi intervals
    print >>text, 'Entries :', np.unique(positive.data)
    positive.data = map(lambda x: mpi(x), positive.data)
    positive = positive.todok()
    sympy_positive = SparseMatrix(size, size, positive)
    L, D = LDL(sympy_positive)
    D = D._smat.values()
    delta = [x.delta/x.mid for x in D]
    print >>text, "The smallest diagonal element in LDL': ", min(D)
    print >>text, 'Ratio largest/smallest :  ', max(D)/min(D)
    print >>text, "Maximal relative delta around diagonal elements: ", max(delta)
    print >>text, '\n'
text.close()
