import math
import itertools
import numpy as np
from igraph import *

def geographicalDiversity(g, p):
    """
    g: Graph
    p: Path between two nodes : p = [L, N], where L =  links id, N = nodes id
    return: The minimun distance between any node members of vector p and the shortest path
    """
    N = p[1]
    s = N[0]
    d = N[len(N) -1]
    shortestPath = g.get_shortest_paths(s, d)[0]
    m = float('inf')
    cond = False
    for node in N:
        for node_s in shortestPath:
            if(g.vertex_disjoint_paths(node, node_s, neighbors = "ignore") != 0): #If exist a path between node and node_s
                distance = g.shortest_paths_dijkstra(node, node_s)[0]
                if(distance < m):
                    m = distance
                    cond = True
    
    if(cond):
        return m
    else:
        raise Exception('There is no path between any node members of vector p and the shortest path')

def getAllSimplePaths(g, s, d, visited, partial = [], result= []):
    """
    g: Graph
    s: Source vertex
    d: Destination vertex
    visited: List of booleans each represent if the vertex is visited, type: numpy array
    partial: Partial path
    result: actual result
    return: A list of all the simple paths between vertex s and d
    """
    partial.append(s)
    visited[s] = True
    if(s == d):
        result.append(partial)
        return result
    neighbors = g.neighbors(s)
    for n in neighbors:
        partial_aux = partial.copy()
        visited_aux = visited.copy()
        if (not visited[n]):
            result = getAllSimplePaths(g, n, d, visited_aux, partial_aux, result)
            visited[n] = True

    return result

def effectiveGeographicalPathDiversity(g, s, d, l):
    """
    g: Graph
    s: Source node
    d: Destination node
    l: lambda value
    return:
    """
    k = 0
    visited = [False] * g.vcount()
    simplePaths = getAllSimplePaths(g, s, d, visited)
    for path in simplePaths:
        k += geographicalDiversity(g, path)
    
    return 1 - math.exp(- l * k)

def totalGraphGeographicalDiversity(g, l):
    """
    g: Graph
    l : lambda value for effective geographical path diversity
    return: The average of effective geographical path diversity of all vertex pairs within g.
    """
    nodes = list(range(g.vcount()))
    pairs = list(itertools.permutations(nodes, 2)) #This list does not contain pairs (v,v)
    sum = 0
    for p in pairs:
        sum += effectiveGeographicalPathDiversity(g, p[0], p[1], l)
    return sum / len(pairs)

def compensatedTotalGeographicalGraphDiversity(g, l):
    """
    g: Graph
    l : lambda value for total graph geographical path diversity
    return: The total graph geographical diversity weighted by the total number of edges of g
    """
    e = g.ecount()
    return totalGraphGeographicalDiversity(g, l) * e

def functionality(g, i, v, attack):
    """
    g: Graph
    i: Attack numer i that eliminates vertex v_i
    v: Vertex
    attack: Attack sequence
    """
    v = g.vcount()
    if(i > v):
        raise Exception('i can not be greater than the total number of vertices')
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

def functionalityLoss(g, v, attack):
    """
    g: Graph
    v: Vertex
    k: Number of attacks
    attack: Attack sequence
    return: The sum of the differences of vertex functionality before and after the attack over the sequence of k attacks
    """
    sum = 0
    k = len(attack)
    for i in range(1, k):
        sum += functionality(g, i -1, v, attack) - functionality(g, i , v, attack)
    return sum

def globalFunctionalityLoss(g, attack):
    """
    g: Graph
    attack: Attack sequence
    edge: Condition that indicates if the attack is for edges
    """
    sum = 0
    v = g.vcount()
    for i in range(v):
        if i not in attack:
            sum += functionalityLoss(g, i, attack)
        
def shortestTemporalDistance(g, s, d, t2):
    """
    g: Graph
    s: Source vertex
    d: Destination vertex
    t2: 
    return:
    """
    if(g.vertex_disjoint_paths(s, d, neighbors = "ignore") != 0): #If there is a way between vertex s and d
        path_len = g.shortest_paths_dijkstra(s, d)[0]
        if(path_len > t2):
            return float('inf')
        else:
            return path_len

def temporalEfficiency(g, t2):
    """
    g: Graph
    t2:
    return:
    """
    sum = 0
    v = g.vcount()
    for i in range(v):
        for j in range(v):
            if(i != j):
                if(g.vertex_disjoint_paths(i, j, neighbors = "ignore") != 0):
                    sum += (1 / shortestTemporalDistance(g, i, j, t2))
    return (1 / (v * (v - 1))) * sum

def deltaEfficiency(g, i):
    """
    g: Graph
    i: Vertex
    return: Efficiency of vertex i
    """
    sum = 0
    v = g.vcount()
    for j in range(v):
        if (j != i):
            if(g.vertex_disjoint_paths(i, j, neighbors = "ignore") != 0):
                    sum += (1 / shortestTemporalDistance(g, i, j, float('inf')))
    return (1 / (v - 1)) * sum

def fragility(g, t2):
    """
    g: Graph
    t2:
    return
    """
    sum = 0
    v = g.vcount()
    for k in range(v):
        g_k = g.copy()
        v_k = list(range(v))
        v_k.remove(k)
        #TODO: remove all links connected to node k
        for i in v_k:
            for j in v_k:
                if(i != j and g.vertex_disjoint_paths(i, j, neighbors = "ignore") != 0 and g_k.vertex_disjoint_paths(i, j, neighbors = "ignore") != 0):
                    sum += (1 / shortestTemporalDistance(g, i, j, t2)) - (1 / shortestTemporalDistance(g_k, i, j, t2))
    return (1 / (v * (v - 1) * (v - 2)) * sum) 

def dynamicFragility(g, t2):
    """
    g: Graph
    t2:
    """
    return 1 - fragility(g, t2)