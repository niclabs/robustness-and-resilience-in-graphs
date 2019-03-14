from igraph import *
import random
import numpy as np
from scipy import linalg

def vertexLoad(g, v, n):
    neighbors = g.neighbors(v)
    nDegree = g.degree(neighbors)
    s = sum(nDegree)
    return (s * g.degree(v)) ** n

def randomWalk(g, s, t, l):
   """
   Random walk between s and t
   l must contain s
   g: Graph
   s: Source vertex
   t: Target
   l: List
   """
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

def entropyRankFromMatrix(m):
    """
    m: Adjancecy matrix
    """
    w, vl, vr = linalg.eig(m, left = True)
    l = len(vl)
    s = np.empty(l)
    for i in range(l):        
        #absvl = np.absolute(vl[i])
        #absvr = np.absolute(vr[i])
        absvl = np.abs(vl[i].real)
        absvr = np.abs(vr[i].real)
        s[i] = np.dot(absvl, absvr)

    s = s/sum(s)
    return s


def entropyRank(g):
    """
    g: Graph
    """
    return entropyRankFromMatrix(g.get_adjacency().data)

def freeEnergyRank(g, e):
    """
    g: Graph
    e: Small number that replaces zeros in the adjacency matrix
    """
    m = np.array(g.get_adjacency().data, dtype = float)
    m[m == 0] = e
    return entropyRankFromMatrix(m)



    


       
#g = Graph([(0,1), (0,2), (2,3), (3,4), (4,2), (2,5), (5,0), (6,3), (5,6)])
g = Graph([(0,1),(0,3), (1,3), (3,3)])
#print(g.get_adjacency())
#print(randomWalkBetweenness(g))
#print(criticality(g, 0))
#print(freeEnergyRank(m))

#Ejemplo matriz paper 24
m = np.array([[0,1,1,1,0,0,0,0], [1,0,1,1,0,0,0,0], [1,1,0,1,1,0,0,0], [1,1,1,0,0,0,0,0], [0,0,0,0,0,0,1,0], [0,0,0,0,1,0,0,0], [0,0,0,0,0,1,0,1], [0,1,0,0,0,0,0,0]])
#print(entropyRankFromMatrix(m))
print(freeEnergyRank(m, 0.003))
