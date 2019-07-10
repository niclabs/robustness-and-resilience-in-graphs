from igraph import *
import random
import numpy as np
import math
from scipy import linalg
import sys
import sympy as sym

def vertexLoad(g, v, n):
    neighbors = g.neighbors(v)
    nDegree = g.degree(neighbors)
    s = sum(nDegree)
    return (s * g.degree(v)) ** n

def randomWalk(g, i, t, s= 0):
   """
   g: Graph
   i: Source vertex
   t: Target
   s: Seed for randomize
   return: Random walk between s and t, list where each element is a vertex
   """
   if (g.vertex_disjoint_paths(i, t, neighbors = "ignore") == 0): #If there is no path between s and t
       return []
   if(s):
       random.seed(s)
   l = [i]
   actual = i
   while(actual != t):
       neighbors = g.neighbors(actual, mode='OUT')
       if len(neighbors) == 0:
           actual = i
           l = [i]
       else:
           n = random.choice(neighbors)
           l.append(n)
           actual = n
   return l

def vertexWalkToEdgesWalk(g, l):
    """
    g: Graph
    l: List, random walk
    return: Converts a list of vertices into a list of edges
    """
    new_l = []
    for i in range(len(l) - 1):
        new_l.append(g.get_eid(l[i], l[i+1]))
    return new_l

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

def criticality(g, v, w, edge = False, s = 0):
    """
    g: Graph
    v: Vertex or edge for criticality
    w: String that represent the weight attributte in the graph
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

def entropyRankFromMatrix(m, i):
    """
    m: Adjancecy matrix
    i: vertex
    return: the entropy rank of vertex i
    """
    w, u, v = linalg.eig(m, left = True)
    w_index = np.argmax(w)
    u_vector = u[w_index]
    v_vector = v[w_index]
    #Normalize u_vector, sum(u_vector_i) = 1
    sum_u = np.sum(u_vector)
    u_vector = u_vector / sum_u
    #Normalize v_vector, sum(u_vector_i * v_vector_i) = 1
    v_vector = v_vector / np.dot(u_vector, v_vector)

    return u_vector[i] * v_vector[i]

def entropyRank(g , i):
    """
    g: Graph
    i: vertex
    return: the entropy rank of vertex i
    """
    m = np.array(g.get_adjacency().data)
    return entropyRankFromMatrix(m, i)

def freeEnergyRank(g, i, e):
    """
    g: Graph
    i: verterx
    e: Small number that replaces zeros in the adjacency matrix
    return: the free energy rank of vertex i
    """
    m = np.array(g.get_adjacency().data)
    m[m == 0] = e
    return entropyRankFromMatrix(m, i)

def bridgeness(g, i, j):
    """
    g: Graph
    i: First component of the edge
    j: Second component of the edge
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

def coveringDegree(g, v):
    """
    g: Graph
    v: Vertex
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


def mcv(m,result, partial,  directed):
    """
    m: Adjacency matrix of a graph
    return: The set of minimal vertex covers of m
    """
    if (len(m) == 0):
        return []
    if (np.count_nonzero(m) == 0):
        partial.sort()    
        if not partial in result:
            result.append(partial)
        return result
    else:
        s = np.sum(m, axis = 0)
        max = np.where(s > 0)[0] #Array that contains all vertex with incident edges
        for u in max:
            mc = np.copy(m)
            partialCopy = list(partial)
            if not directed:
                #For each neighbor of u, delete incident edges
                nu = np.where(mc[u] == 1)[0]
                for i in nu:
                    mc[:,i] = 0
            mc[:,u] = 0
            partialCopy.append(u)
            result = mcv(mc, result, partialCopy, directed)
        return result

def MCV(g):
    """
    g: Graph
    return: The set of minimum vertex covers of G
    """
    m = g.get_adjacency().data
    covers = mcv(m, [], [], g.is_directed())
    result = []
    min = len(covers[0])
    for cover in covers:
        l = len(cover)
        if l < min:
            min = l
    for cover in covers:
        if(len(cover) == min):
            result.append(cover)
    return result

def coveringIndex(g, v):
    """
    g: Graph
    v: Vertex
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

def sensitivity(g, f, w, t, k):
    """
    g:Graph
    f: Centrality function: f(M), where M is the adjacency matrix
    w: Name of the weight attribute in the graph g
    t:
    k:
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

def weightedMatrix(g, w):
    """
    g: Graph
    w: Strings that represent the name of the edge weight attribute in the graph
    return: The weighted adjacency matrix of graph g
    """
    m = g.get_adjacency().data
    n = len(m)
    for i in range(n):
        for  j in range(n):
            if(m[i][j] != 0):
                e_id = g.get_eid(i,j)
                e_weight = g.es[e_id].attributes()[w]
                m[i][j] = e_weight
    return m
