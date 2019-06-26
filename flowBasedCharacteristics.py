import numpy as np
import math
from igraph import *
import scipy

def electricalNodalRobustness(g, i, attribute):
    """
    g: Graph
    i: vertex
    attribute: Name of the edge attribute that contains the flow
    return: The electrical nodal robustness of a vertex i
    """
    edges = g.adjacent(i) #id's of the edges vertex i is incident on
    L = len(edges)
    f = np.zeros(L)
    sum = 0   
    for j in range(L):
        fi = g.es[edges[j]].attributes()[attribute]
        f[j] = fi
    
    sumf = np.sum(f)
    k = 0
    for e in edges:
        edge = g.get_edgelist[e:e+1][0]
        incidentVertex = edge[1]
        Li = g.degree(incidentVertex, mode="OUT")
        pi = (f[k] / sumf)
        sum += ( 1 / (L * Li)) * pi * math.log(pi)
        
    return - sum

def relativeAreaIndex(g, v, w, f, u0, ut):
    """
    g: Graph
    v: vertex
    w: Weight function over a parameter u
    f: Maximum flow at v given u, use: f(v, u)
    u0:
    ut:
    return: The relative area index of vertex v
    """
    funnum = lambda u : w(u) * (f(v, u0) - f(u))
    funden = lambda u : w(u) * (f(v, u0))
    num = 0
    den = 0
    for x in range(u0, ut + 1):
        num += scipy.misc.derivative(funnum, x)
        den += scipy.misc.derivative(funden, x)
    return num / den