from igraph import *
import numpy as np
from auxiliaryFunctions import *

def percentageOfNoncriticalNodes(t, g, g_nodes, c_nodes, a, b):
    """
    t: Tolerance threshold
    g: Power grid    
    g_nodes: Generation nodes
    c_nodes: Consumer nodes
    a: Tolerance parameter of nodes
    b: Tolerance parameter of links
    """
    pass

def robustnessValue(g, rank=None):
    """
    g: Graph
    rank: Rank function, it gives a list of nodes in order of importance, default: Pagerank
    """
    rank_list = []
    n = g.vcount()
    if rank is None:
        rank_list = np.array(g.pagerank())
        rank_list = rank_list.argsort()[::-1][:n]
    else:
        rank_list = rank(g)
    
    g_aux = g.copy()
    #Remove nodes and calculate the size of giant connected component
    for i in range(n - 1):
        vertex = rank_list[i]
        neighbors = g_aux.neighbors(vertex, mode='OUT')
        #Delete edges from vertex
        for target in neighbors:
            g_aux.delete_edges([(vertex, target)])

        #Get size of giant component
        comp = g_aux.components()
        rank_list[i] = max(comp.sizes())
    rank_list[n-1] = 0

    return np.sum(rank_list) / n

def treeness(g):
    """
    g: Graph
    """
    Gc = getDAG(g)
    Gc_r = getDAG(g, reverse=True)
    hf = H_f(Gc)
    hb = H_f(Gc_r)
    return hf- hb / max(hf, hb)

def entropy(graph, gamma= 0.5, x_pos=None, y_pos=None):
    """
    gamma: Weight of the node criticality, weight of edge criticality 1 - gamma
    x_pos: Parameter that indicates the name attribute for x position of vertices
    y_pos: Parameter that indicates the name attribute for y position of vertices
    If x_pos and y_pos are not indicated will be set randomly and returned with the final result
    """
    if x_pos is None or y_pos is None:
        x_pos = np.random.rand(graph.vcount())
        y_pos = np.random.rand(graph.vcount())    
    else:
        x_pos = graph.vs[x_pos]
        y_pos = graph.vs[y_pos]
    #Calculate node criticality
    CV = criticality(graph, graph.vcount())
   
    edge_correlation_factor = np.zeros((graph.vcount(), graph.vcount()))
    #Calculate edge correlation factor matrix
    for i in range(graph.vcount()):
        for j in range(graph.vcount()):
            if i == j:
                edge_correlation_factor[i][j] = 1
            else:
                edge_correlation_factor[i][j] = max(CV[i], CV[j]) / (CV[i] + CV[j])
    
    adjacency = graph.get_adjacency().data
    #Change adjacency matrix
    for i in range(graph.vcount()):
        for j in range(g.vcount()):
            if i == j:
                adjacency[i][j] = 1
    
    average_transmission_efficiency = np.zeros(graph.vcount())
    #Calculate average transmission efficiency vector
    for i in range(graph.vcount()):
        acc = 0
        for j in range(graph.vcount()):
            if i!= j:
                distance = 0 #TODO: Calculate
                acc += 1 / distance
        average_transmission_efficiency[i] = (2 * acc) / (graph.vcount() * (graph.vcount() - 1))
    
    #Average transmission efficiency matrix
    I = np.zeros((graph.vcount(), graph.vcount()))
    for i in range(graph.vcount()):
        for j in range(graph.vcount()):
            I[i][j] = average_transmission_efficiency[i]
    
    #Calculate edge criticality matrix
    CE = I * adjacency * edge_correlation_factor

    #Calculate comprehensive critical degree of nodes
    CS = np.zeros(graph.vcount())
    for i in range(graph.vcount()):
        acc = 0
        i_neighbors = graph.neighbors(i)
        for j in i_neighbors:
            acc += CE[i][j]
        CS[i] = (gamma * CV[i]) + (1 - gamma) * acc / len(i_neighbors)

    #Calculate critical coefficient of nodes
    S = CS / np.sum(CS)
    result = - np.sum( S * np.log(S))
    return result, x_pos, y_pos

def globalResilienceAnalysis():
    pass

def robustnessMeasureR(graph, ranking_function=None):
    """
    ranking_function: Function that returns a list of ranked nodes in graph, use: ranking_function(graph), if ranking_function is None, then the nodes are ranked by degree
    """
    if ranking_function is None:
        degrees = graph.degree()
        ranking = np.flip(np.argsort(degrees))
    else:
        ranking = ranking_function(graph)
    acc = 0
    n = graph.vcount()
    for vertex in ranking:
        graph.delete_vertices(vertex)
        s = sizeMaxComponent(graph)
        acc += s
    try:
        result = acc / n
    except ZeroDivisionError:
        result = None
    return result

def globalConnectivity(graph):
    n = graph.vcount()
    try:
        result = sizeMaxComponent(graph) / n
    except ZeroDivisionError:
        result = None
    return result

def localConnectivity(graph):
    giant_component = getGiantComponent(graph)
    try:
        result = localConnectivityAux(graph) / localConnectivityAux(giant_component)
    except ZeroDivisionError:
        result = None
    return result
