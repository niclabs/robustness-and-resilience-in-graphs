import math
import itertools
import numpy as np
from igraph import *
from auxiliaryFunctions import *

def geographicalDiversity(g, p=[0,1]):
    """
    g: Graph
    p: Vertex path between two nodes
    return: The minimun distance between any node members of vector p and the shortest path
    """
    s = p[0]
    d = p[len(p) -1]
    shortestPath = g.get_shortest_paths(s, d)[0]
    m = float('inf')
    cond = False
    for node in p:
        for node_s in shortestPath:
            if(node != node_s and g.vertex_disjoint_paths(node, node_s, neighbors = "ignore") != 0): #If exist a path between node and node_s
                distance = g.shortest_paths_dijkstra(node, node_s)[0][0]
                if(distance < m):
                    m = distance
                    cond = True
    
    if(cond):
        return m
    else:
        raise Exception('There is no path between any node members of vector p and the shortest path')

def effectiveGeographicalPathDiversity(g, s= 0, d=1, l= 1):
    """
    g: Graph
    s: Source node, default = 0
    d: Destination node, default = 1
    l: Constant to weight the utility of alternatie paths, default = 1
    """
    k = 0
    visited = [False] * g.vcount()
    print(f"Calculando {s}, {d}")
    simplePaths = getAllSimplePaths(g, s, d, visited)
    print("Listo")
    for path in simplePaths:
        k += geographicalDiversity(g, path)
    
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

def deltaEfficiency(g, i=0):
    """
    g: Graph
    i: Vertex, default = 0
    return: Efficiency of vertex i
    """
    sum = 0
    v = g.vcount()
    for j in range(v):
        if (j != i):
            if(g.vertex_disjoint_paths(i, j, neighbors = "ignore") != 0):
                    sum += (1 / shortestTemporalDistance(g, i, j, float('inf')))
    return (1 / (v - 1)) * sum

def fragility(g, t2=1):
    """
    g: Graph
    t2: Limit of temporal distance, default= 1
    return: The average delta efficiency over single-vertex removals
    """
    sum = 0
    v = g.vcount()
    for k in range(v):
        g_k = g.copy()
        v_k = list(range(v))
        v_k.remove(k)
        #Remove all links connected to node k
        edges = get_edges(g_k, k)
        for edge in sorted(edges, reverse= True):
                g_k.delete_edges(edge)
        for i in v_k:
            for j in v_k:
                if(i != j and g.vertex_disjoint_paths(i, j, neighbors = "ignore") != 0 and g_k.vertex_disjoint_paths(i, j, neighbors = "ignore") != 0):
                    sum += (1 / shortestTemporalDistance(g, i, j, t2)) - (1 / shortestTemporalDistance(g_k, i, j, t2))
    return (1 / (v * (v - 1) * (v - 2)) * sum) 

def dynamicFragility(g, t2=1):
    """
    g: Graph
    t2: Limit of temporal distance, default= 1
    return: The average delta efficiency over all vertices, weighted by the failure probability of each vertex
    """
    return 1 - fragility(g, t2)