from igraph import *
import numpy as np

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

def definition523(g, k, degreeProduct= True):
    """
    g: Graph
    k: Number of edge removal
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

def maxDegreeProduct(g):
    """
    g: Graph
    return: Edge id of the edge with the max degree product
    """
    m = g.ecount()
    degree = np.zeros(m)
    for i in range(m):
        s = g.get_edgelist()[i: i+1][0][0] #Source vertex of the edge
        d = g.get_edgelist()[i: i+1][0][1] #Destination vertex of the edge
        degree[i] = g.degree(s) * g.degree(d)

    return degree.index(max(degree))

def maxEdgeBetweenness(g):
    """
    g: Graph
    return: Edge id of the edge with the max edge betweenness
    """
    bet = g.edge_betweenness(directed=g.is_directed())
    return bet.index(max(bet))