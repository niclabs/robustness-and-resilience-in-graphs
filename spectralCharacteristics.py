import numpy as np
import math
#from scipy import heaviside
import scipy.sparse.linalg as sp
from igraph import *
from auxiliaryFunctions import *

def reconstructabilityCoefficient(g):
    """
    g: Graph
    """
    A = np.array(g.get_adjacency().data)
    A = A.astype(float)
    eigenvalues, eigenvectors = np.linalg.eigh(A) #Eigen values are not necessarily ordered, eigenvectors[:,i] is the eigenvector corresponding to the eigenvalues[i]
    eigenvectors = normalize(eigenvectors)
    result = 0
    A_i = eigenvectors * np.diag(eigenvalues) * np.transpose(eigenvectors)
    A_i = np.heaviside(A_i, 0.5)
    for j in range(len(eigenvalues)):
        result += 1
        eigenvalues[j] = 0
        new_a = eigenvectors * np.diag(eigenvalues) * np.transpose(eigenvectors)
        #new_a -= 0.5
        new_a = np.heaviside(new_a, 0.5)
        if not np.array_equal(new_a, A_i): #Compare
            return result
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
    eigenvalues = np.abs(eigenvalues)
    sum = 0
    for i in range(k):
        sum += (eigenvectors[i][v] ** 2) + math.sinh(eigenvalues[i])
    return sum


def generalizedRobustnessIndex(g, k= 1):
    """
    g: Graph
    k: Number of eigenvalues to include into the approximation, default = 1
    return:
    """
    A = np.array(g.get_adjacency().data)
    A = A.astype('float64')
    w, v = sp.eigsh(A, k, which='LM' )
    SC = np.matmul(v**2, np.sinh(w))
    v[v <= 0] = None
    SC[SC <= 0] = None

    log_v = np.log10(v[:, 0])
    log_v = np.nan_to_num(log_v)
    log_SC = 0.5 * np.log10(SC)
    log_SC = np.nan_to_num(log_SC)
    sinh_w = np.sinh(w[0])
    log_w = 0
    if(sinh_w > 0):
        log_w = 0.5 * np.log10(np.sinh(w[0]))

    d = log_v - log_SC + log_w
    r_k = (np.sum(d ** 2) / len(A)) ** 0.5
    return r_k

def localNaturalConnectivity(g): # auxiliary measure
    """
    g: Graph
    """
    A = np.array(g.get_adjacency().data)
    eigenvalues = np.linalg.eig(A)[0]
    n = len(eigenvalues)
    return np.log(closedWalkNumber(g)/n)

def closedWalkNumber(g): # auxiliary
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
    comp = g.components(mode='WEAK')
    sum = 0
    for c in comp:
        aux_m = m.copy()
        aux_m = changeMatrix(m, c)
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
