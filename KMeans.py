#!/usr/bin/env python
# kmeans.py using any of the 20-odd metrics in scipy.spatial.distance
# kmeans_sample 2 pass, first sample sqrt(N)

from __future__ import division
import numpy as np
from scipy.spatial.distance import cdist
from scipy.sparse import issparse


def kmeans(X, centres, delta=.001, maxiter=10, metric="euclidean", p=2, verbose=1):

    """
    :param X: (N, dim) matrix. Can be sparse.
    :param centres: (k, dim) initial centres. E.g. random.sample(X, k).
    :param delta: Relative error. Iterate until the avg distance to centres in within delta.
    :param metric: any of the 20-odd in scipy.spatial.distance.
    :param p: for minkowski metric -- local mod cdist for 0 < p < 1 too.
    
    :return:
        centres, (k, dim)
        Xtocentre: each X -> its nearest centre, ints N -> k.
        distances, N
    """

    if not issparse(X):
        X = np.asanyarray(X)
    centres = centres.todense() if issparse(centres) \
        else centres.copy()
    N, dim = X.shape
    k, cdim = centres.shape

    if dim != cdim:
        raise ValueError("kmeans: X %s and centres %s must have the same number of columns" % (
            X.shape, centres.shape))
    if verbose:
        print "kmeans: X %s  centres %s  delta=%.2g  maxiter=%d  metric=%s" % (
            X.shape, centres.shape, delta, maxiter, metric)

    allx = np.arange(N)
    prev_dist = 0
    for jiter in range(1, maxiter + 1):

        D = cdist_sparse(X, centres, metric=metric, p=p)  # |X| x |centres|
        xtoc = D.argmin(axis=1)  # X -> nearest centre
        distances = D[allx, xtoc]
        avdist = distances.mean()  # median ?

        if verbose >= 2:
            print "kmeans: av |X - nearest centre| = %.4g" % avdist
        if (1 - delta) * prev_dist <= avdist <= prev_dist \
                or jiter == maxiter:
            break

        prev_dist = avdist
        for jc in range(k):  # (1 pass in C)
            c = np.where(xtoc == jc)[0]
            if len(c) > 0:
                centres[jc] = X[c].mean(axis=0)
    if verbose:
        print "kmeans: %d iterations  cluster sizes:" % jiter, np.bincount(xtoc)
    if verbose >= 2:
        r50 = np.zeros(k)
        r90 = np.zeros(k)
        for j in range(k):
            dist = distances[xtoc == j]
            if len(dist) > 0:
                r50[j], r90[j] = np.percentile(dist, (50, 90))
        print "kmeans: cluster 50 % radius", r50.astype(int)
        print "kmeans: cluster 90 % radius", r90.astype(int)
        # scale L1 / dim, L2 / sqrt(dim) ?
    return centres, xtoc, distances


def kmeans_sample(X, k, nsample=0, **kwargs):
    """ 2-pass kmeans, fast for large N:
        1) kmeans a random sample of nsample ~ sqrt(N) from X
        2) full kmeans, starting from those centres
    """
    # merge w kmeans ? mttiw
    # v large N: sample N^1/2, N^1/2 of that
    # seed like sklearn ?
    N, dim = X.shape
    if nsample == 0:
        nsample = max(2 * np.sqrt(N), 10 * k)
    Xsample = random_sample(X, int(nsample))
    pass1centres = random_sample(X, int(k))
    sample_centres = kmeans(Xsample, pass1centres, **kwargs)[0]
    return kmeans(X, sample_centres, **kwargs)


def cdist_sparse(X, Y, **kwargs):
    """ -> |X| x |Y| cdist array, any cdist metric
        X or Y may be sparse -- best csr
    """
    # todense row at a time, v slow if both v sparse
    sxy = 2 * issparse(X) + issparse(Y)
    if sxy == 0:
        return cdist(X, Y, **kwargs)
    d = np.empty((X.shape[0], Y.shape[0]), np.float64)
    if sxy == 2:
        for j, x in enumerate(X):
            d[j] = cdist(x.todense(), Y, **kwargs)[0]
    elif sxy == 1:
        for k, y in enumerate(Y):
            d[:, k] = cdist(X, y.todense(), **kwargs)[0]
    else:
        for j, x in enumerate(X):
            for k, y in enumerate(Y):
                d[j, k] = cdist(x.todense(), y.todense(), **kwargs)[0]
    return d


def random_sample(X, n):
    """ random.sample of the rows of X
        X may be sparse -- best csr
    """
    sampleix = random.sample(xrange(X.shape[0]), int(n))
    return X[sampleix]


def nearestcentres(X, centres, metric="euclidean", p=2):
    """ each X -> nearest centre, any metric
            euclidean2 (~ withinss) is more sensitive to outliers,
            cityblock (manhattan, L1) less sensitive
    """
    D = cdist(X, centres, metric=metric, p=p)  # |X| x |centres|
    return D.argmin(axis=1)


def Lqmetric(x, y=None, q=.5):
    # yes a metric, may increase weight of near matches; see ...
    return (np.abs(x - y) ** q).mean() if y is not None \
        else (np.abs(x) ** q).mean()


class Kmeans:
    """ km = Kmeans( X, k= or centres=, ... )
        in: either initial centres= for kmeans
            or k= [nsample=] for kmeans_sample
        out: km.centres, km.Xtocentre, km.distances
        iterator:
            for jcentre, J in km:
                clustercentre = centres[jcentre]
                J indexes e.g. X[J], classes[J]
    """

    def __init__(self, X, k=0, centres=None, nsample=0, **kwargs):
        self.X = X
        if centres is None:
            self.centres, self.Xtocentre, self.distances = kmeans_sample(
                X, k=k, nsample=nsample, **kwargs)
        else:
            self.centres, self.Xtocentre, self.distances = kmeans(
                X, centres, **kwargs)

    def __iter__(self):
        for jc in range(len(self.centres)):
            yield jc, (self.Xtocentre == jc)


if __name__ == "__main__":
    import random
    import sys
    from time import time

    N = 10000
    dim = 10
    ncluster = 10
    kmsample = 100  # 0: random centres, > 0: kmeans_sample
    kmdelta = .001
    kmiter = 10
    metric = "cityblock"  # "chebyshev" = max, "cityblock" L1,  Lqmetric
    seed = 1

    exec ("\n".join(sys.argv[1:]))  # run this.py N= ...
    np.set_printoptions(1, threshold=200, edgeitems=5, suppress=True)
    np.random.seed(seed)
    random.seed(seed)

    print "N %d  dim %d  ncluster %d  kmsample %d  metric %s" % (
        N, dim, ncluster, kmsample, metric)
    X = np.random.exponential(size=(N, dim))
    # cf scikits-learn datasets/
    t0 = time()
    if kmsample > 0:
        centres, xtoc, dist = kmeans_sample(X, ncluster, nsample=kmsample,
                                           delta=kmdelta, maxiter=kmiter, metric=metric, verbose=2)
    else:
        random_centres = random_sample(X, ncluster)
        centres, xtoc, dist = kmeans(X, random_centres,
                                     delta=kmdelta, maxiter=kmiter, metric=metric, verbose=2)
    print "%.0f msec" % ((time() - t0) * 1000)
