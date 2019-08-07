from igraph import *
import random
import numpy as np
import math
from scipy import linalg
import sys
import sympy as sym
from auxiliaryFunctions import *

def vertexLoad(g, v=0, n=1):
    """
    g: Graph
    v: Vertex, default=0
    n: Power, default = 1
    return: The degree of vertex v multiplied by the sum of the degrees of its neighbors elevated
            to a power
    """
    neighbors = g.neighbors(v)
    nDegree = g.degree(neighbors)
    s = sum(nDegree)
    return (s * g.degree(v)) ** n

def randomWalkBetweenness(g, edge = False, seed = 0):
    """
    g: Graph
    edge: Boolean that indicates if we count edges
    seed: Seed for randomize
    return: An array 'a' where a[i] is the randomWalkBetweenness of element i, the element can be a vertex or an edge
    """
    if edge:
        sum = np.zeros(g.ecount(), dtype = int)
    else:
        sum = np.zeros(g.vcount(), dtype = int)
    for s in range(g.vcount()):
        for t in range(g.vcount()):
            if s != t and g.vertex_disjoint_paths(s,t, neighbors = "ignore") != 0: #Si hay un camino
                    walk = randomWalk(g, s, t, seed)
                    if edge: #Si contamos arcos
                        walk = vertexWalkToEdgesWalk(g, walk)
                        for e in range(g.ecount()):
                            sum[e] += walk.count(e)
                    else:
                        for v in range(g.vcount()): #Contar la aparición de cada vértice
                            sum[v] += walk.count(v)
    return sum

def criticality(g, v=0, w='weight', edge = False, s = 0):
    """
    g: Graph
    v: Vertex or edge for criticality, default=0
    w: String that represent the weight attributte in the graph, default='weight'
    edge: Condition that indicates if v is an edge
    s: Seed for randomize
    return: The random-walk betweenness of an element normalized by the weight of the element
    """
    betweenness = randomWalkBetweenness(g, edge, seed=s) #Array
    if edge:
        weight = g.es[v].attributes()[w]
    else:
        weight = g.vs[v].attributes()[w]
    if weight == 0:
        return 0
    return betweenness[v] / weight

def entropyRank(g , i=0):
    """
    g: Graph
    i: vertex, defualt=0
    return: the entropy rank of vertex i
    """
    m = np.array(g.get_adjacency().data)
    return entropyRankFromMatrix(m, i)

def freeEnergyRank(g, i=0, e=0.01):
    """
    g: Graph
    i: verterx, default = 0
    e: Small number that replaces zeros in the adjacency matrix, default = 0.01
    return: the free energy rank of vertex i
    """
    m = np.array(g.get_adjacency().data)
    m[m == 0] = e
    return entropyRankFromMatrix(m, i)

def bridgeness(g, i=0, j=1):
    """
    g: Graph, assumes exists edge (0,1) by default
    i: First component of the edge, default=0
    j: Second component of the edge default=1
    """
    if(not g.are_connected(i,j)):
        raise Exception("Edge doesn't exist")
    cliques = g.cliques() #List of tuples, each one is a clique
    Si = 1
    Sj = 1
    Se = 1
    for tuples in cliques:
        l = np.array(tuples)
        if(i in l):
            if(len(l) > Si):
                Si = len(l)
            if (j in l):
                if(len(l) > Se):
                    Se  = len(l)
        if(j in l):
            if(len(l) > Sj):
                Sj = len(l)

    return math.sqrt(Si * Sj) / Se

def coveringDegree(g, v=0):
    """
    g: Graph
    v: Vertex, defualt = 0
    return: The number of minimal vertex cover that contains v
    """
    m = np.array(g.get_adjacency().data)
    if(np.sum(m, axis=0)[v] == 0): #v has no incident edges
        return 0
    covers = mcv(m, [], [], g.is_directed())
    result = 0
    for cover in covers:
        if v in cover:
            result += 1
    return result

def coveringIndex(g, v=0):
    """
    g: Graph
    v: Vertex, default=0
    """
    m = g.get_adjacency().data
    c = len(mcv(m, [], [], g.is_directed()))
    mVertCov = MCV(g)
    a = 0
    for cover in mVertCov:
        if v in cover:
            a += 1
    b = coveringDegree(g, v)
    return a + b/c

def sensitivity(g, f, t, k, w='weight'):
    """
    g:Graph
    f: Centrality function: f(M), where M is the adjacency matrix
    t: Number of the argument in the adjacency matrix function, ex: A(t_0, ..., t_n)
    k: Variation matrix is with respect to node k
    w: Name of the weight attribute in the graph g, default='weight'
    return: The sensitivity matrix
    """
    h = sym.Symbol('h')
    M = sym.Symbol('M')
    A = weightedMatrix(g, w)
    n = len(A)
    dAdt = np.zeros([n, n])
    deg = g.degree(t)
    for i in range(n):
        for j in range(n):
            Aij = A[i][j]
            if( i == k or j == k):               
                value = sym.limit((Aij * (1 + (h / deg)) - Aij) / h, h, 0)
            else:
                value = sym.limit((Aij- Aij) / h, h, 0)
            dAdt[i][j] = value
    dfdA = sym.diff(f, M)
    return dfdA * dAdt