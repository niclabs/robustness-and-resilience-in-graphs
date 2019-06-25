import numpy as np
import math
from igraph import *

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