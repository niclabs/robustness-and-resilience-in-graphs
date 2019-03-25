from igraph import *
import random
import numpy as np
import math
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
    w, u, v = linalg.eig(m, left = True)
    l = len(u)
    s = np.empty(l)
    for i in range(l):
        absu = np.absolute(u[:,i])
        absv = np.absolute(v[:,i])
        #absu = np.abs(u[:,i].real)
        #absv = np.abs(v[:,i].real)
        s[i] = np.dot(absu, absv)

    s = s/sum(s)
    return s


def entropyRank(g):
    """
    g: Graph
    """
    m = np.array(g.get_adjacency().data)
    return entropyRankFromMatrix(m)

def freeEnergyRank(g, e):
    """
    g: Graph
    e: Small number that replaces zeros in the adjacency matrix
    """
    m = np.array(g.get_adjacency().data)
    m[m == 0] = e
    return entropyRankFromMatrix(m)

def bridgeness(g, i, j):
    """
    g: Graph
    i: First component of the edge
    j: Second component of the edge
    """
    if(not g.are_connected(i,j)):
        raise Exception("Edge ("+ str(i) + ", " + str(j)+ ") doesn't exist")
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

def searchEdge(m, u = 0, vertex= False, directed= False):
    """
    Auxiliary function for coveringDegree
    m: Adjacency matrix
    u: 
    vertex: Condition that indicates if u is given
    return: convenient edge [u v]
    """
    if (np.count_nonzero(m) == 0): #m is zero
        return
    if not vertex:
        #Find u
        sumIncEdges = np.sum(m, axis= 0) #Amount of incident edges of each vertex
        u = np.argmax(sumNeighbors)
    aux_m = m #Auxiliary matrix to verify conditions
    aux_m[:, u] = 0
    if not directed:
        aux_m[u] = 0
    if (np.count_nonzero(m) == 0): #vertex cover is u
        return [u]
         
    #Find v
    uNeighbors = np.where(m[u] == 1)[0]
    v = uNeighbors[0]
    max = 0
    sumColumns = np.sum(m, axis = 0)
    for n in uNeighbors:
        n_i = sumColumns[n]
        if(max < n_i):
            v = n
            max = n_i
    return [u, v]

def coveringDegree(g, v):
    """
    g: Graph
    v: Vertex
    return: The number of minimal vertex cover that contains v
    """   
    result = []
    m = np.array(g.get_adjacency().data)
    if(np.sum(m, axis= 0)[v] == 0): #v has no incident edges
        return 0
    else:
        edge = searchEdge(m, v, vertex=True)
    #While de adjacency matrix isn't zero  
    while(np.count_nonzero(m) != 0):
        #TODO: incluir cuando edge es solo un vertice
        result.append(edge[0])
        result.append(edge[1])
        #Remove all edges incident on u or v
        m[:, edge[0]] = 0
        m[:, edge[1]] = 0
        if not g.is_directed():
            m[edge[0]] = 0
            m[edge[1]] = 0
        edge = searchEdge(m)
        print(edge)
    return result
    #return len(result) #TODO:en lugar de contar acá, tener una variable que vaya guardando la cantidad

def coveringDegree2(g, v):
    """
    g: Graph
    v: Vertex
    return: The number of minimal vertex cover that contains v
    """
    m = np.array(g.get_adjacency().data)
    if(np.sum(m, axis=0)[v] == 0): #v has no incident edges
        return 0
    result = []
    result.append(v)
    if not g.is_directed():
        #TODO: Para todos los vecinos de v borrar todos los arcos incidentes
        pass
    m[:,0] = 0 #Delete edges incident to v
    while (np.count_nonzero(m) != 0):
        u = np.argmax(np.sum(m, axis = 0))
        result.append(u)
        if not g.is_directed():
            #TODO: Para todos los vecinos de u borrar todos los arcos incidentes
            pass
        m[:,u] = 0
    return result






    


       
#g = Graph([(0,1), (0,2), (2,3), (3,4), (4,2), (2,5), (5,0), (6,3), (5,6)])
#g = Graph([(0,1),(0,3), (1,3), (3,3)])
#print(randomWalkBetweenness(g))
#print(criticality(g, 0))

#Ejemplo matriz paper 24
#g = Graph([(0,1), (0,2), (0,3), (1,0), (1,2), (1,3), (2,0), (2,1), (2,3), (2,4), (3,0), (3,1), (3,2), (4,6),(5,4), (6,5), (6,7), (7,1)], directed = True)
#print(entropyRank(g))

g = Graph([(0,1), (1,2) , (1,3), (3,4)])
m = np.array(g.get_adjacency().data)
#print(searchEdge(m,0, vertex=True))
