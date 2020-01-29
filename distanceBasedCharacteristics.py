import math
import itertools
import numpy as np
from igraph import *
import networkx as nx
from auxiliaryFunctions import *

def geographicalDiversity(g, p=[0,1],distance=False):
    """
    g: Graph
    p: Vertex path between two nodes
    distance: Name of the edge attribute that indicates distance
    return: The minimun distance between any node members of vector p and the shortest path, return inf when if there is not a path
    """
    if not distance:
        distance="distance"
        g = generateWeight(g,name=distance)

    s = p[0]
    d = p[-1]

    #Get shortest path
    shortestPath = g.get_shortest_paths(s, d, weights=distance)[0]
    
    m = float('inf')
    for node in p:
        for node_s in shortestPath:
            if(node != node_s and g.vertex_disjoint_paths(node, node_s, neighbors = "ignore") != 0): #If exist a path between node and node_s
                distance = g.shortest_paths_dijkstra(node, node_s)[0][0]
                if(distance < m):
                    m = distance
    if(m == float('inf')):
        return 0
    return m

def effectiveGeographicalPathDiversity(g, s= 0, d=1, l= 1):
    """
    g: Graph
    s: Source node, default = 0
    d: Destination node, default = 1
    l: Constant to weight the utility of alternatie paths, default = 1
    """
    k = 0

    #Get all simple paths
    edge_list = g.get_edgelist()
    n_g = nx.Graph(edge_list)
    simplePaths = list(nx.all_simple_paths(n_g, s,d))

    for path in simplePaths:
        k += geographicalDiversity(g, list(path))
    
    return 1 - math.exp(- l * k)

def totalGraphGeographicalDiversity(g, l=1):
    """
    g: Graph
    l : Constant for effective geographical path diversity, to weight the utility of alternatie paths, default = 1
    return: The average of effective geographical path diversity of all vertex pairs within g.
    """
    nodes = list(range(g.vcount()))
    pairs = list(itertools.permutations(nodes, 2)) #This list does not contain pairs (v,v)
    sum = 0
    for p in pairs:
        sum += effectiveGeographicalPathDiversity(g, p[0], p[1], l)
    return sum / len(pairs)

def compensatedTotalGeographicalGraphDiversity(g, l=1):
    """
    g: Graph
    l : Constant for total graph geographical path diversity, to weight the utility of alternatie paths, default = 1
    return: The total graph geographical diversity weighted by the total number of edges of g
    """
    e = g.ecount()
    return totalGraphGeographicalDiversity(g, l) * e

def functionality(g, attack=[0], v = 0, i = 0):
    """
    g: Graph
    attack: Sequence of vertices that represents the attack, this list must contain at least one vertex
    v: Vertex, default = 0
    i: Attack number i that eliminates vertex v_i, default = 0
    return: Functionality of vertex v after the attack number i that eliminates vertex vi 
    """
    if(len(attack) == 0):
        return -1
    result = np.zeros(i + 1)
    result[0] = 1
    k = 1
    aux = g.copy()
    aux.delete_vertices([attack[0]])
    while(k <= i):
        distance = aux.shortest_paths_dijkstra(v, k)[0]
        degree = aux.degree(k)
        result[k] = result[k - 1] - (1 / ((distance ** 2) * degree)) * result[k - 1]
        k += 1
        aux.delete_edges([attack[k]])
    return result[i]

def functionalityLoss(g, attack=[0], v=0):
    """
    g: Graph
    attack: Attack sequence, this list must contain at least one vertex
    v: Vertex, dafault = 0
    return: The sum of the differences of vertex functionality before and after the attack over the sequence of k attacks
    """
    sum = 0
    k = len(attack)
    for i in range(1, k):
        sum += functionality(g, i -1, v, attack) - functionality(g, i , v, attack)
    return sum

def globalFunctionalityLoss(g, attack=[0]):
    """
    g: Graph
    attack: Attack sequence, this list must contain at least one vertex
    return: The sum of functionality loss over all vertices in the graph that aren't removed in the attack
    """
    sum = 0
    v = g.vcount()
    for i in range(v):
        if i not in attack:
            sum += functionalityLoss(g, attack, i)
    return sum

def temporalEfficiency(g, t2=1):
    """
    g: Graph
    t2: Limit of temporal distance, default= 1
    return: The way in which the efficiency of g degrades after an elimination
    """
    sum = 0
    v = g.vcount()
    for i in range(v):
        for j in range(v):
            if(i != j):
                if(g.vertex_disjoint_paths(i, j, neighbors = "ignore") != 0):
                    sum += (1 / shortestTemporalDistance(g, i, j, t2))
    return (1 / (v * (v - 1))) * sum

def deltaEfficiency(graph, k=None):
    """
    k: Vertex
    return: Delta effiency of graph with respect to node k
    """
    n = graph.vcount()
    if n == 0:
        return None
    if k is None:
        k = random.randint(0, n - 1)
    acc = 0
    v_k = list(range(n))
    v_k.remove(k)
    graph_k = removeLinks(graph, k)
    for i in v_k:
        for j in v_k:
            if i != j:
                distance_g = float(graph.shortest_paths_dijkstra(i,j)[0][0])
                distance_g_k = float(graph_k.shortest_paths_dijkstra(i, j)[0][0])
                acc += (1/distance_g) - (1/distance_g_k)
    try:
        return acc / ((n-1) * (n-2))
    except ZeroDivisionError:
        return None

def fragility(g):
    """
    g: Graph
    return: The average delta efficiency over single-vertex removals
    """
    acc = float(0)
    for k in range(g.vcount()):
        acc += deltaEfficiency(g, k)
    try:
        return acc / g.vcount()
    except ZeroDivisionError:
        return None

def dynamicFragility(g):
    """
    g: Graph
    return: The average delta efficiency over all vertices, weighted by the failure probability of each vertex
    """
    return 1 - fragility(g)