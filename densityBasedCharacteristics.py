from igraph import *
import numpy as np
from auxiliaryFunctions import *

def hubDensity(g):
    """
    g: Graph, each component of the graph is a hub
    return: The average density of subgraphs induced by hubs
    """
    components = g.components()
    n_components = len(components)
    sum = 0
    for c in components:
        sum += len(c)
    return sum / n_components

def definition523(g, k= 1, degreeProduct= True):
    """
    g: Graph
    k: Number of edge removal, default = 1
    degreeProduct: Degree product is the strategy for the removal, if false, edge betweenness
    return: The average fraction of edges in the largest connected component after k edge removals
    """
    aux = g.copy()
    sum = 0
    for q in range(k):
        if degreeProduct:
            id_edge = maxDegreeProduct(aux)
        else:
            id_edge = maxEdgeBetweenness(aux)
        aux.delete_edges(id_edge)
        comp = aux.components()
        max_len = 0
        for c in comp:
            if len(c) > max_len:
                max_len = c
        sum += max_len
    return sum / k