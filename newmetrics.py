from igraph import *
import numpy as np
import math
from auxiliaryFunctions import *

def robustnessIndex(g, rank=None):
    """
    g: Graph
    rank: Rank function, it gives a list of nodes in order of importance, default: Pagerank
    """
    rank_list = []
    n = g.vcount()
    if rank is None:
        rank_list = np.array(g.pagerank())
        rank_list = rank_list.argsort()[::-1][:n]
    else:
        rank_list = rank(g)
    
    g_aux = g.copy()
    # Remove nodes and calculate the size of giant connected component
    for i in range(n - 1):
        vertex = rank_list[i]
        neighbors = g_aux.neighbors(vertex, mode='OUT')
        # Delete edges from vertex
        for target in neighbors:
            g_aux.delete_edges([(vertex, target)])
        # Get size of giant component
        comp = g_aux.components()
        rank_list[i] = max(comp.sizes())
    rank_list[n-1] = 0

    return np.sum(rank_list) / n

def robustnessMeasureR(graph, ranking_function=None):
    """
    ranking_function: Function that returns a list of ranked nodes in graph, use: ranking_function(graph), if ranking_function is None, then the nodes are ranked by degree
    """
    if ranking_function is None:
        degrees = graph.degree()
        ranking = np.flip(np.argsort(degrees))
    else:
        ranking = ranking_function(graph)
    acc = 0
    n = graph.vcount()
    for vertex in ranking:
        # Delete vertex
        graph.delete_vertices(vertex)
        # Calculate the size of giant connected component
        comp = graph.components()
        acc += max(comp.sizes())
    try:
        result = acc / n
    except ZeroDivisionError:
        result = None
    return result

def RCB(g):
    """
    """
    # Check connected graph
    comp = g.components()
    if len(comp) > 1:
        return None
    m = g.ecount()
    n = g.vcount()
    cospanning_trees = cospanningTrees(g)
    n_t = len(cospanning_trees)
    acc = 0
    for i in range(m):
        n_t_e = 0
        for cosp_tree in cospanning_trees:
            if i in cosp_tree:
                n_t_e += 1
        CB_e_i = n_t_e / ((n - 1) * n_t)
        if CB_e_i != 0:
            acc += CB_e_i * math.log2(CB_e_i)
    return -acc / math.log2(m)