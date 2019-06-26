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
    return: The edge betweenness of edge (s,d)
    """
    sum = 0
    edgeId = g.get_eid(s, d)
    v = g.vcount()
    for i in range(v):
        for j in range(v):
            if(i != j and g.vertex_disjoint_paths(i, j, neighbors = "ignore") != 0):
                randWalk = randomWalk(g, i, j, [i])
                edgeList = vertexWalkToEdgesWalk(g, randWalk)
                sum += edgeList.count(edgeId)
    return sum

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