import numpy as np
import math
from scipy import heaviside
from igraph import *
from auxiliaryFunctions import *

def reconstructabilityCoefficient(g):
    """
    g: Graph
    """
    A = np.array(g.get_adjacency().data)
    eigenvalues, eigenvectors = np.linalg.eig(A) #Eigen values are not necessarily ordered, eigenvectors[:,i] is the eigenvector corresponding to the eigenvalues[i]
    eigenvalues, eigenvectors = sortEigenValuesVectors(eigenvalues, eigenvectors)
    result = 0

    for j in range(1, len(eigenvalues)):
        eigenvalues[j] = 0
        new_a = eigenvectors * np.diag(eigenvalues) * np.transpose(eigenvectors)
        new_a = heaviside(new_a, 0.5)
        if(not np.array_equal(new_a, A)): #Compare
            return result
        result += 1
    return result

def normalizedSubgraphCentrality(g, v= 0, k=1):
    """
    g: Graph
    v: vertex, default = 0
    k: Number of eigenvalues to include into the approximation, default = 1
    return: the normalized subgraph centrality of vertex v
    """
    A = np.array(g.get_adjacency().data)
    eigenvalues, eigenvectors = np.linalg.eig(A)
    eigenvalues, eigenvectors = sortEigenValuesVectors(eigenvalues, eigenvectors)
    sum = 0
    for i in range(k):
        sum += (eigenvectors[i][v] ** 2) + math.sinh(eigenvalues[i])
    return sum

def generalizedRobustnessIndex(g, k=1):
    """
    g: Graph
    k: Number of eigenvalues to include into the approximation, default = 1
    return:
    """
    A = np.array(g.get_adjacency().data)
    eigenvalues, eigenvectors = np.linalg.eig(A)
    eigenvalues, eigenvectors = sortEigenValuesVectors(eigenvalues, eigenvectors)
    sum = 0
    v = g.vcount()
    for i in range(v):
        sum += (math.log(eigenvectors[1][i]) - (math.log( math.sinh(eigenvalues[1]) ** -0.5) + math.log(normalizedSubgraphCentrality(g, i, k)) / 2 ))  ** 2
    return math.sqrt (sum / v)

def localNaturalConnectivity(g):
    """
    g: Graph
    """
    A = np.array(g.get_adjacency().data)
    eigenvalues = np.linalg.eig(A)[0]
    n = len(eigenvalues)
    return np.log(closedWalkNumber(g)/n)

def closedWalkNumber(g):
    """
    g: Graph
    return: A weighted sum of the number of closed walks in a graph
    """
    A = np.array(g.get_adjacency().data)
    eigenvalues = np.linalg.eig(A)[0]
    exp = np.exp(eigenvalues)
    return np.sum(exp)

def redundancyOfAlternativePaths(g):
    """
    g: Graph
    return:
    """
    return closedWalkNumber(g)

def naturalConnectivity(g):
    """
    g: Graph
    """
    m = g.get_adjacency().data
    eigenvalues = np.linalg.eig(m)[0]
    n = len(eigenvalues)
    comp = g.components()
    sum = 0
    for c in comp:
        aux_m = m.copy()
        aux_m[np.ix_(c,c)]
        l_i = naturalConAux(aux_m)
        n_i = len(c)
        sum += n_i * np.exp(l_i)
    
    return np.log(sum / n)

def subgraphCentrality(g):
    """
    g: Graph
    return: The sum of closed walks for all vertices
    """
    return np.log(closedWalkNumber(g))

def normalizedLocalNaturalConnectivity(g):
    """
    g: Graph   
    """
    return localNaturalConnectivity(g)