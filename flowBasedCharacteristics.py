import numpy as np
import math
from igraph import *
from scipy.misc import derivative
from auxiliaryFunctions import *

def electricalNodalRobustness(g, i=0, attribute='flow'):
    """
    g: Graph
    i: vertex, default = 0
    attribute: Name of the edge attribute that contains the flow, default= 'flow'
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

def relativeAreaIndex(g, w=weightFunction, f=maxFlow, v=0, u0= 0, ut= 1):
    """
    g: Graph
    w: Weight function over a parameter u
    f: function that represents the maximum flow at v given u, use: f(v, u)
    v: vertex, default = 0
    u0: First limit of integral, default = 0
    ut: Second limit of the integral, default = 1
    return: The relative area index of vertex v
    """
    funnum = lambda u : w(u) * (f(v, u0) - f(v,u))
    funden = lambda u : w(u) * (f(v, u0))
    num = 0
    den = 0
    for x in range(u0, ut + 1):
        num += derivative(funnum, x)
        den += derivative(funden, x)
    return num / den