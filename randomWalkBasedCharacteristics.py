from igraph import *
from auxiliaryFunctions import *

def networkCriticality(g, w=False, vertices=True, edges=False):
    """
    g: Graph
    w: String, name of the atribute for link weight, default = False, if false, weights are set randomly
    vertices: Sum criticality of vertices
    edges: Sum criticality of edges
    return: The sum of criticalities over all elements of G (vertices, edges or both)
    """
    if not w:
        g = generateWeight(g, edge=True, vertex=True)
        w = 'weight'
    sum = 0
    v = g.vcount()
    for i in range(v):
        if vertices:
            sum += criticalityOfVertex(g, i, w)
        if edges:
            for j in range(v):
                if (i != j and g.get_eid(i, j, directed=False, error=False) != -1): #If there is a link
                    sum += criticalityOfEdge(g, i, j, w)
    return sum