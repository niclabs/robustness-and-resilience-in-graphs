from igraph import *

def definition100(g, k):
    """
    g: Graph
    k: Number of edge removal
    return: The average fraction of edges in the largest connected component after 
            k edge removals
    """
    m = g.vcount()
    aux = g.copy()
    sum = 0
    for q in range(m):
        



    