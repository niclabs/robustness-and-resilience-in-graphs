from igraph import *
from auxiliaryFunctions import *

def networkCriticality(g, w='weight', vertices=True, edges=False):
    """
    g: Graph
    w: String, name of the atribute for link weight, default = 'weight'
    vertices: Sum criticality of vertices
    edges: Sum criticality of edgfes
    return: The sum of criticalities over all elements of G (vertices, edges or both)
    """
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