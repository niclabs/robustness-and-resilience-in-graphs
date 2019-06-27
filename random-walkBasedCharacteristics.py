from igraph import *
from characteristicsForRankingElements import randomWalk, vertexWalkToEdgesWalk

def criticalityOfEdge(g,i, j, w):
    """
    g: Graph
    i: Source vertex of the edge
    j: Incident vertex
    w: String, name of the atribute for link weight
    return: Criticality of edge e, it assumes that edge exists
    """
    link_id = g.get_eid(i,j)
    edgeBet = g.edge_betweenness(directed=g.is_directed())[link_id]
    return edgeBet / g.es[link_id].attributes()[w]

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

def networkCriticality(g, w, vertices=True, edges=False):
    """
    g: Graph
    w: String, name of the atribute for link weight
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