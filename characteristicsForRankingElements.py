from igraph import *
import random
import numpy as np
import math
from scipy import linalg
import sys

def vertexLoad(g, v, n):
    neighbors = g.neighbors(v)
    nDegree = g.degree(neighbors)
    s = sum(nDegree)
    return (s * g.degree(v)) ** n

def randomWalk(g, s, t, l, seed = 0):
   """
   Random walk between s and t
   l must contain s
   g: Graph
   s: Source vertex
   t: Target
   l: List
   seed: Seed for randomize
   """
   if(seed):
       random.seed(seed)
   neighbors = g.neighbors(s)
   if len(neighbors) == 0:
       return l
   else:
       n = random.choice(neighbors)
       if n == t:
           l.append(n)
           return l
       else:
           l.append(n)
           return randomWalk(g, n, t, l)

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

def randomWalkBetweenness(g, edge = False):
    if edge:
        sum = np.zeros(g.ecount(), dtype = int)
    else:
        sum = np.zeros(g.vcount(), dtype = int)
    for s in range(g.vcount()):
        for t in range(g.vcount()):
            if s != t and g.vertex_disjoint_paths(s,t, neighbors = "ignore") != 0: #Si hay un camino
                    walk = randomWalk(g, s, t, [s])
                    while walk[len(walk) - 1] != t: #Para que el camino llegue a t (sea valido)
                        walk = randomWalk(g, s, t, [s])
                    #Luego de encontrado el camino
                    if edge: #Si contamos arcos
                        walk = vertexWalkToEdgesWalk(g, walk)
                        for e in range(g.ecount()):
                            sum[e] += walk.count(e)
                    else:
                        for v in range(g.vcount()): #Contar la aparición de cada vértice
                            sum[v] += walk.count(v)

    return sum

def criticality(g,v, edge = False):
    """
    g: Graph
    v: Vertex or edge for criticality
    edge: Condition that indicates if v is an edge
    """
    betweenness = randomWalkBetweenness(g, edge)
    s = sum(betweenness)
    if s == 0:
        return 0
    return betweenness[v]/s

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

def sensitivity(g, f, v, w):
    """
    g:Graph
    f: Centrality function: f(d_1, ..., d_n) with n the number of vertex
    v: Vertex
    w: Vertex
    return: Sensitivity of node v with respect to the node w
    """
    return
