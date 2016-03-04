#!/usr/bin/env sage -python
"""
Find rational LU decomposition for our matrices.

Using sparse sage matrices over rationals (QQ), and built in LU decomposition.

Each matrix takes 6-12 hours to decompose. Sorted diagonal entries saved in
LU_rational folder in rational and floating point form. All diagonal entries
for all matrices are positive, hence matrices are positive-definite.
"""

import numpy as np
from sage.all import matrix, QQ, Rational, Set, parallel
from datetime import datetime


def get_matrix(filename):
    """
    Import matrix from file.

    Entries assumed to be integers.
    """
    f = open(filename)
    v = np.array(eval(f.readline()), dtype=int)
    f.close()
    size = int(max(v[:, 0])) + 1
    dct = dict([((x[0], x[1]), Rational(x[2])) for x in v])
    return matrix(QQ, size, size, dct, sparse=True)


@parallel(ncpus=4)
def process_matrix(number):
    """ Parallel processing worker. """
    start = datetime.now()
    positive = get_matrix(
        'results/CHOLMOD_permuted/permuted{}.txt'.format(number))
    _, L = positive.LU(pivot='nonzero', format='compact')
    diag = sorted(list(Set(L.diagonal())))
    f = open('results/LU_rational/diagonal{}.txt'.format(number), 'w')
    for d in diag:
        print >>f, d, '\n', float(d), '\n'
    f.close()
    text = open('results/lu_sage.txt', 'a')
    print >>text, 'Domain : ', number
    print >>text, '{}. Smallest diagonal entry: '.format(number), float(diag[0])
    print >>text, '{}. Time: '.format(number), datetime.now() - start
    print >>text, '\n'
    text.close()

text = open('results/lu_sage.txt', 'w')
text.close()
#list(process_matrix(range(1, 21)))
