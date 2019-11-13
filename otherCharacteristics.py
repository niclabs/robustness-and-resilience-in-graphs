from igraph import *
import numpy as np
from auxiliaryFunctions import *


def effectiveGraphResistance(g, weight=False):
    """
    g: Graph
    w: Name of the attribute that contains the weight of each edge, default = 'weight'
    return: The sum of all effective resistances between all pairs in a network
    """
    if not weight:
        g = generateWeight(g) #Edge weight
        q = g.laplacian(weights = 'weight')
    else:
        q = g.laplacian(weights = weight)
    q_plus = np.linalg.pinv(q)
    v = g.vcount()
    acc = 0
    for i in range(v):
        for j in range(v):
            acc += q_plus[i][i] - 2 * q_plus[i][j] + q_plus[j][j]
    return acc

def viralConductance(g):
    """
    g: Graph
    return: the average fraction of infected nodes
    """
    m = g.get_adjacency().data
    eigenvalues = np.linalg.eig(m)[0]
    max_e = np.amax(eigenvalues)
    max_e = int(np.real(max_e))
    acc = 0

    for s in range(max_e + 1):
        acc += y(g, s) 
    mean = acc / (max_e + 1)  
    return max_e * mean