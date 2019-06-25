from igraph import *

def criticalityOfEdge(g,i, j, w):
    """
    g: Graph
    i: Source vertex of the edge
    j: Incident vertex
    w: String, name of the atribute for link weight
    return: Criticality of edge e
    """
    link_id = g.get_eid(i,j)
    return edgeBetweenness(g, i, j) / g.es[link_id].attributes()[w]

def criticalityOfVertex(g, v, w):
    """
    g: Graph
    v: Vertex
    w: String, name of the atribute for link weight
    return: Criticality of vertex v
    """
    sum = 0
    neighbors = g.neighbors(v, mode= 'IN')
    for n in neighbors:
        link_id = g.get_eid(n, v)
        sum += g.es[link_id].attributes()[w]
    return g.betweenness(v, directed=g.is_directed(), weights=w) / sum

def edgeBetweenness(g, s, d):
    """
    g: Graph
    s: Source vertex
    d: Destination vertex
    return:
    """
    #TODO
    return 0

def networkCriticality(g, w, onlyEdges=False, both=False):
    """
    g: Graph
    w: String, name of the atribute for link weight
    onlyEdges: If false, then only vertices
    both
    return: The sum of criticalities over all elements of G (vertices, edges or both)
    """
    sum = 0
    v = g.vcount()
    for i in range(v):
        if (not onlyEdges):
            sum += criticalityOfVertex(g, i, w)
        if both:
            for j in range(v):
                if (i != j):
                    sum += criticalityOfEdge(g, i, j, w)
    return sum