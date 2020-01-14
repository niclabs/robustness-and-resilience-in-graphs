from igraph import *
import numpy as np
from auxiliaryFunctions import *

def percentageOfNoncriticalNodes(t, g, g_nodes, c_nodes, a, b):
    """
    t: Tolerance threshold
    g: Power grid    
    g_nodes: Generation nodes
    c_nodes: Consumer nodes
    a: Tolerance parameter of nodes
    b: Tolerance parameter of links
    """
    pass

def robustnessValue(g, rank=None):
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
    #Remove nodes and calculate the size of giant connected component
    for i in range(n - 1):
        vertex = rank_list[i]
        neighbors = g_aux.neighbors(vertex, mode='OUT')
        #Delete edges from vertex
        for target in neighbors:
            g_aux.delete_edges([(vertex, target)])

        #Get size of giant component
        comp = g_aux.components()
        rank_list[i] = max(comp.sizes())
    rank_list[n-1] = 0

    return np.sum(rank_list) / n

def treeness(g):
    """
    g: Graph
    """
    Gc = getDAG(g)
    Gc_r = getDAG(g, reverse=True)
    hf = H_f(Gc)
    hb = H_f(Gc_r)
    return hf- hb / max(hf, hb)


