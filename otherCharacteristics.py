from igraph import *
import numpy as np


def effectiveGraphResistance(g, weight):
    """
    g: Graph
    w: Name of the attribute that contains the weight of each edge
    return: The sum of all effective resistances between all pairs in a network
    """
    q = g.laplacian(weights = weight)
    q_plus = np.linalg.pinv(q)
    v = g.vcount()
    sum = 0
    for i in range(v):
        for j in range(v):
            sum += q_plus[i][i] - 2 * q_plus[i][j] + q_plus[j][j]
    return sum



