# Statistical Analysis of Network Data Project, ENSAE 2015.
# Robin Vogel, Clement Puppo, Alain Soltani.
# ---------------------

import networkx as nx
import scipy as sp
from scipy.sparse import spdiags, csr_matrix


# 1. DER Algorithm.
# ---------------------

def average_partition_degree(subgraph, nodes=None, weight=None):

    """Average weighted degree within partition."

    Parameters
    ----------
    subgraph : graph of partition to compute average degree on.
        A NetworkX partition

    weight : matrix of weights, accounting for the average of
        the empirical measures of sequences x that start at i.
        scipy.sparse Matrix
    """

    weighted_sum, normalization = csr_matrix((1, weight.shape[1])), sum(subgraph.degree.values())

    for n, deg in subgraph.degree.items():

        weighted_sum += subgraph.degree[n].value() * weight.getrow(n)

    return weighted_sum / normalization


def der_algorithm(g, l, k, nodelist=None, weight='weight'):

    """
    Parameters
    ----------
    g : graph
       A NetworkX graph

    l : Walk length.

    k : number of components.

    nodelist : list, optional
       The rows and columns are ordered according to the nodes in nodelist.
       If nodelist is None, then the ordering is produced by G.nodes().

    weight : string or None, optional (default='weight')
       The edge data key used to provide each value in the matrix.
       If None, then each edge has weight 1.
    """

    # a: Adjacency matrix.

    if nodelist is None:
        nodelist = g.nodes()
    a = nx.to_scipy_sparse_matrix(g, nodelist=nodelist, weight=weight,
                                  dtype='float')
    n, m = a.shape

    # Di : Inverse degree matrix.

    di = spdiags(1.0/sp.array(a.sum(axis=1).flat), [0], n, n)

    # T : lists of transition matrices, from 1 to l.

    t = [(di * a) ** (i + 1) for i in range(l)]

    # w_i : distribution corresponding to the average of
    # the empirical measures of sequences x that start at i.

    w = sum(t) / l

    # Random graph partition.

    # TODO: Select p, q accordingly to the p,q-SBM analytic bounds.
    p, q = 0.25, 0.01

    rpg = nx.random_partition_graph(
        [g.number_of_nodes() / k for i in range(k)], p, q)

    print rpg

    """# DER Algorithm.

    while "THE SETS P DO NOT CHANGE":

        # 1. For all partitions, construct each mean and store it.

        mean_list = [average_partition_degree(sub_rpg, w) for sub_rpg in rpg.graph['partition']]

        # 2. Update all partitions accordingly to the maximisation step.

        for sub_rpg in rpg.graph['partition']:
            # sub_rpg = ...

    return rpg.graph['partition']"""

if __name__ == "__main__":

    G = nx.karate_club_graph()
    print der_algorithm(G, 10, 2)




