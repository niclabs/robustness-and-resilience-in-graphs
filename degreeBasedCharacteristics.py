import math
import numpy as np
from igraph import *

def degreeEntropy(g):
    """
    g: Graph
    """
    sum = 0
    for i in range(1, g.vcount()):
        p_i = getProbabilityDegree(g, i)
        sum += p_i * math.log(p_i)
    return -sum

def relativeEntropy(g):
    """
    g: Graph
    """
    n = g.vcount()
    pk = getDegreeDistribution(g)
    sum = 0
    for i in range(n):
        sum += pk[i] * math.log(n * pk[i])
    return sum

def getProbabilityDegree(g, k, myMode='ALL'):
    """
    g: Graph
    k: 
    myMode: degree distribution mode, it can be 'ALL', 'IN' or 'OUT'
    return: The probability p(k) of the degree distribution
    """
    h = g.degree_distribution(mode= myMode) #Degree distribution
    bins = list(h.bins())
    acc = 0
    for b in bins:
        min = math.floor(b[0])
        if(min <= k):
            acc += b[2]
        else:
            break   
    return acc / g.vcount()

def getDegreeDistribution(g, myMode='ALL'):
    """
    g: Graph
    myMode: degree distribution mode, it can be 'ALL', 'IN' or 'OUT'
    return: List with the degree distribution
    """
    h = g.degree_distribution(mode= myMode) #Degree distribution
    bins = list(h.bins())
    p_k = np.zeros(len(bins))
    i = 0
    for b in bins:
        p_k[i] += b[2]
        i += 1   
    return p_k / g.vcount()