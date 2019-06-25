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

def viralConductance(g):
    """
    g: Graph
    return: the average fraction of infected nodes
    """
    m = g.get_adjacency().data
    eigenvalues = np.linalg.eig(m)[0]
    max = np.amax(eigenvalues)
    sum = 0

    for i in range(max + 1):
        sum += y(g, i) 
    mean = sum / (max + 1)  
    return max * mean

def probis(g, i, s):
    """
    g: Graph
    i: Vertex
    s: State
    return: probability that node i is infected at steady state s.
    """
    neighbors = g.neighbors(i)
    sum = 0
    for j in neighbors:
        sum += probis(g, j, s)
    return sum / (s + sum)

def y(g, s):
    """
    g: Graph
    s:
    return: Fraction of infected nodes
    """
    sum = 0
    v = g.vcount()
    for i in range(v):
        sum += probis(g,i, s)
    return 1 / (v * sum)
