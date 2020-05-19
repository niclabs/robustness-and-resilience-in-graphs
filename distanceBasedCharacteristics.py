import math
import itertools
import numpy as np
from igraph import *
import networkx as nx
from auxiliaryFunctions import *

def geographicalDiversity(g, p=None, x_position=None, y_position=None):
    """
    g: Graph
    p: Vertex path between two nodes
    x_position: Name of the vertex attribute that indicates x position of vertex
    y_position: Name of the vertex attribute that indicates y position of vertex
    return: The minimun distance between any node members of vector p and the shortest path, return inf when if there is not a path
    """
    n = g.vcount()
    if n < 2:
        return None
    if x_position is None:
        g = generatePositions(g, 'x', 'y')
        x_position = 'x'
        y_position = 'y'
    distance = 'distance'
    g = generateDistance(g, x_position, y_position, distance)  
    if p is None:
        vertices = list(range(n))
        s = random.choice(vertices)
        vertices.remove(s)
        d = random.choice(vertices)
        p = randomWalk(g, s, d)
    s = p[0]
    d = p[-1]
    #Get shortest path
    shortestPath = g.get_shortest_paths(s, d, weights=distance)[0]  
    m = float('inf')
    for node in p:
        for node_s in shortestPath:
            if node != node_s and node not in shortestPath:
                x1 = g.vs[x_position][node]
                y1 = g.vs[y_position][node]
                x2 = g.vs[x_position][node_s]
                y2 = g.vs[y_position][node_s]
                distance_i = math.sqrt(((x1- x2) ** 2) + (y1- y2) ** 2)
                if(distance_i < m):
                    m = distance_i
    if(m == float('inf')):
        return 0
    return m

def effectiveGeographicalPathDiversity(g, s= 0, d=1, l= 1, x_position=None, y_position=None):
    """
    g: Graph
    s: Source node, default = 0
    d: Destination node, default = 1
    l: Constant to weight the utility of alternatie paths, default = 1
    x_position: Name of the vertex attribute that indicates x position of vertex
    y_position: Name of the vertex attribute that indicates y position of vertex
    """
    if x_position is None:
        g = generatePositions(g, 'x', 'y')
        x_position = 'x'
        y_position = 'y'
    k = 0
    #Get all simple paths
    edge_list = g.get_edgelist()
    n_g = nx.Graph(edge_list)
    simplePaths = list(nx.all_simple_paths(n_g, s,d))

    for path in simplePaths:
        k += geographicalDiversity(g, list(path), x_position, y_position)   
    return 1 - math.exp(- l * k)

def totalGraphGeographicalDiversity(g, l=1, x_position=None, y_position=None):
    """
    g: Graph
    l : Constant for effective geographical path diversity, to weight the utility of alternatie paths, default = 1
    x_position: Name of the vertex attribute that indicates x position of vertex
    y_position: Name of the vertex attribute that indicates y position of vertex
    return: The average of effective geographical path diversity of all vertex pairs within g.
    """
    if x_position is None:
        g = generatePositions(g, 'x', 'y')
        x_position = 'x'
        y_position = 'y'
    nodes = list(range(g.vcount()))
    pairs = list(itertools.permutations(nodes, 2)) #This list does not contain pairs (v,v)
    sum = 0
    for p in pairs:
        sum += effectiveGeographicalPathDiversity(g, p[0], p[1], l, x_position=x_position, y_position=y_position)
    return sum / len(pairs)

def compensatedTotalGeographicalGraphDiversity(g, l=1, x_position=None, y_position=None):
    """
    g: Graph
    l : Constant for total graph geographical path diversity, to weight the utility of alternatie paths, default = 1
    x_position: Name of the vertex attribute that indicates x position of vertex
    y_position: Name of the vertex attribute that indicates y position of vertex
    return: The total graph geographical diversity weighted by the total number of edges of g
    """
    if x_position is None:
        g = generatePositions(g, 'x', 'y')
        x_position = 'x'
        y_position = 'y'
    e = g.ecount()
    return totalGraphGeographicalDiversity(g, l, x_position=x_position, y_position=y_position) * e

def functionality(g, attack=None, v = None, l = None): # auxiliary measure
    """
    g: Graph
    attack: Sequence of vertices that represents the attack, this list must contain at least one vertex
    v: Vertex, default = random
    m: Attack number m that eliminates vertex v_i, default = random
    return: Functionality of vertex v after the attack number i that eliminates vertex vi 
    """
    if g.vcount() < 2:
        return None
    vertices = list(range(g.vcount()))
    #Set random default values
    if v is None:
        v = random.choice(vertices)
    if attack is None:
        total_attack_length = random.randint(1, g.vcount() - 1)
        vertices.remove(v)
        attack = random.sample(vertices, total_attack_length)
    if l is None:
        l = random.randint(0, len(attack))
    #Compute measure
    result = np.zeros(l + 1)
    result[0] = float(1)
    k = 1
    aux = g.copy()
    while(k <= l):
        v_i = attack[k - 1]
        distance = aux.shortest_paths_dijkstra(v_i, v)[0][0]
        degree = aux.degree(v)
        if degree == 0:
            result[k] = 0
        else:
            result[k] = result[k - 1] - (result[k - 1] / ((distance ** 2) * degree)) 
        aux = removeLinks(aux, v_i)
        k += 1
    return result[l]

def functionalityLoss(g, attack=None, v=None, m=None): # auxiliary measure
    """
    g: Graph
    attack: Attack sequence, this list must contain at least one vertex
    v: Vertex, dafault = 0
    m: Amount of attacks to consider, not required
    return: The sum of the differences of vertex functionality before and after the attack over the sequence of k attacks
    """
    if g.vcount() < 2:
        return None
    vertices = list(range(g.vcount()))
    #Set random default values
    if v is None:
        v = random.choice(vertices)
    if attack is None:
        total_attack_length = random.randint(1, g.vcount() - 1)
        vertices.remove(v)
        attack = random.sample(vertices, total_attack_length)
    if m is None:
        m = len(attack)
    f_0 = functionality(g, attack, v, 0)
    f_m = functionality(g, attack, v, m)
    return f_0 - f_m

def globalFunctionalityLoss(g, attack=None):
    """
    g: Graph
    attack: Attack sequence, this list must contain at least one vertex
    return: The sum of functionality loss over all vertices in the graph that aren't removed in the attack
    """
    n = g.vcount()
    if n < 2:
        return None
    vertices = list(range(g.vcount()))
    #Set random default values
    if attack is None:
        total_attack_length = random.randint(1, g.vcount())
        attack = random.sample(vertices, total_attack_length)   
    acc = 0
    n = g.vcount()
    for vertex in range(n):
        if vertex not in attack:
            m = len(attack)
        else:
            m = attack.index(vertex) - 1
        if m < 0:
            pass
        else:
            acc += functionalityLoss(g, attack, vertex, m)
    return acc

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

def fragility(g): # auxiliary
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

def vulnerability(graph, attack_function=attack_edges):
    """
    attack_function: Function that damages graph by deleting vertices, edges or both, returns a graph, use: attack_function(graph)
    """
    performance_normal = performance(graph)
    attacked_graph = attack_function(graph)
    performance_damage = performance(attacked_graph)
    try:
        result = (performance_normal - performance_damage) / performance_damage
    except ZeroDivisionError:
        result = None
    return result
