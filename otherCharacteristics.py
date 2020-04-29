from igraph import *
import numpy as np
from auxiliaryFunctions import *

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
    # Calculate node criticality
    CV = criticality(graph, graph.vcount())
   
    edge_correlation_factor = np.zeros((graph.vcount(), graph.vcount()))
    # Calculate edge correlation factor matrix
    for i in range(graph.vcount()):
        for j in range(graph.vcount()):
            if i == j:
                edge_correlation_factor[i][j] = 1
            else:
                edge_correlation_factor[i][j] = max(CV[i], CV[j]) / (CV[i] + CV[j])
    
    adjacency = graph.get_adjacency().data
    # Change adjacency matrix
    for i in range(graph.vcount()):
        for j in range(graph.vcount()):
            if i == j:
                adjacency[i][j] = 1
    
    average_transmission_efficiency = np.zeros(graph.vcount())
    # Calculate average transmission efficiency vector
    for i in range(graph.vcount()):
        acc = 0
        for j in range(graph.vcount()):
            if i!= j:
                distance = 0 #TODO: Calculate
                acc += 1 / distance
        average_transmission_efficiency[i] = (2 * acc) / (graph.vcount() * (graph.vcount() - 1))
    
    # Average transmission efficiency matrix
    I = np.zeros((graph.vcount(), graph.vcount()))
    for i in range(graph.vcount()):
        for j in range(graph.vcount()):
            I[i][j] = average_transmission_efficiency[i]
    
    # Calculate edge criticality matrix
    CE = I * adjacency * edge_correlation_factor

    # Calculate comprehensive critical degree of nodes
    CS = np.zeros(graph.vcount())
    for i in range(graph.vcount()):
        acc = 0
        i_neighbors = graph.neighbors(i)
        for j in i_neighbors:
            acc += CE[i][j]
        CS[i] = (gamma * CV[i]) + (1 - gamma) * acc / len(i_neighbors)

    # Calculate critical coefficient of nodes
    S = CS / np.sum(CS)
    result = - np.sum( S * np.log(S))
    return result, x_pos, y_pos

def effectiveGraphResistance(g, weight=False):
    """
    g: Graph
    w: Name of the attribute that contains the weight of each edge, default = 'weight'
    return: The sum of all effective resistances between all pairs in a network
    """
    if not weight:
        g = generateWeight(g) #Edge weight
        q = g.laplacian(weights = 'weight')
    else:
        q = g.laplacian(weights = weight)
    q_plus = np.linalg.pinv(q)
    v = g.vcount()
    acc = 0
    for i in range(v):
        for j in range(v):
            acc += q_plus[i][i] - 2 * q_plus[i][j] + q_plus[j][j]
    return acc

def viralConductance(g):
    """
    g: Graph
    return: the average fraction of infected nodes
    """
    m = g.get_adjacency().data
    eigenvalues = np.linalg.eig(m)[0]
    max_e = np.amax(eigenvalues)
    max_e = int(np.real(max_e))
    acc = 0

    for s in range(max_e + 1):
        acc += y(g, s) 
    mean = acc / (max_e + 1)  
    return max_e * mean

def RCB(g):
    """
    """
    # Check connected graph
    comp = g.components()
    if len(comp) > 1:
        return None
    m = g.ecount()
    n = g.vcount()
    cospanning_trees = cospanningTrees(g)
    n_t = len(cospanning_trees)
    acc = 0
    for i in range(m):
        n_t_e = 0
        for cosp_tree in cospanning_trees:
            if i in cosp_tree:
                n_t_e += 1
        CB_e_i = n_t_e / ((n - 1) * n_t)
        if CB_e_i != 0:
            acc += CB_e_i * math.log2(CB_e_i)
    return -acc / math.log2(m)