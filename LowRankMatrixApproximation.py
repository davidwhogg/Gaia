'''
This file is part of the Gaia project.
Copyright 2012 David W. Hogg
'''

import numpy as np
import scipy.sparse as sp

def makeFakeNormalMatrix(npar):
    '''
    Build a npar x npar covariance matrix by random methods.
    Guaranteed to produce non-singular matrices.

    Bugs: Various hard-coded little things; the method is magic.
    '''
    C = np.zeros((npar, npar))
    for i in range(npar+2):
        v = np.random.normal(size=npar)
        C += np.outer(v, v)
    return C

def makeFakeGaiaNormalMatrix(nstars=10):
    '''
    Build a sparse realistic matrix with the same structure as the
    expected Gaia normal matrix, although with many fewer stars and
    attitude knots!

    Bugs: Doesn't use sparse at all yet; I am writing this on a plane
    with no docs.  Many magic numbers.
    '''
    nknots = nstars / 10
    knotsperstar = 70
    npar = 5
    N = np.matrix(np.zeros((nstars * npar + nknots, nstars * npar + nknots)))
    rho = 0.1
    i = nstars * npar
    for knot in range(nknots):
        N[i,i] += 1
        if i > 0:
            N[i-1,i] += 0.5
        if i+1 < nknots:
            N[i+1,i] += 0.5
        i += 1
    i = 0
    for star in range(nstars):
        N[i:i+npar,i:i+npar] += makeFakeNormalMatrix(npar)
        knots = nstars * npar + np.random.randint(nknots, size=(knotsperstar))
        for j in range(npar):
            for k in knots:
                val = rho * N[i, i] * N[k,k]
                N[k, i] += val
                N[i, k] += val
            i += 1
    return N

def main():
    '''
    Step through the project steps
    '''
    M = makeFakeGaiaNormalMatrix()
    print M
    exactC = la.inverse(M)
    return None

if __name__ == '__main__':
    import cProfile as cp
    cp.run('main()')

