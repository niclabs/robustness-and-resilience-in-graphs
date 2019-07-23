from igraph import *

def criticalityOfEdge(g, i= 0, j= 1, w='weight'):
    """
    g: Graph
    i: Source vertex of the edge, default = 0
    j: Incident vertex, default = 1
    w: String, name of the atribute for link weight, default = 'weight'
    return: Criticality of edge e, it assumes that edge exists, default edge (0, 1)
    """
    link_id = g.get_eid(i,j)
    edgeBet = g.edge_betweenness(directed=g.is_directed())[link_id]
    return edgeBet / g.es[link_id].attributes()[w]

def criticalityOfVertex(g, v= 0, w= 'weight'):
    """
    g: Graph
    v: Vertex, default = 0
    w: String, name of the atribute for link weight, default = 'weight'
    return: Criticality of vertex v
    """
    sum = 0
    neighbors = g.neighbors(v, mode= 'IN')
    for n in neighbors:
        link_id = g.get_eid(n, v)
        sum += g.es[link_id].attributes()[w]
    return g.betweenness(v, directed=g.is_directed(), weights=w) / sum

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