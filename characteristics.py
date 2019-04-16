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
        max = np.where(s > 0)[0] #Array that contains al vertex with incident edges
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

def splittingNumber(g, k):
    """
    g: Graph
    k: 
    return: The average number of edges that need to be removed to break  g into k connected components
    """

    numOfLoops = 0 #TODO: k??
    dev = 1 #Se puede cambiar
    givenComponents = 0 #Int

    i = 0
    numOfComponents = 0
    cond = True
    while (cond):
        #TODO: make all links functioning
        numOfLinks[i] = 0 #TODO: definirlo

        while (numOfComponents < givenComponents):
            #TODO:choose randomly an uninterrupted link
            numOfLinks[i] += 1
            #TODO: interrup the link
            #TODO:compute numOfComponents

            if(i > 1 and i > numOfLoops):
                mean1 = np.mean(numOfLinks[1:i+1]) #mean value of numOfLinnks[1:i]
                mean2 = np.mean(numOfLinks[1:i]) #mean value of numOfLinnks[1:i - 1]

                if(abs(mean1 - mean2) < dev):
                    cond = False
        i += 1
    return i



    


       
#g = Graph([(0,1), (0,2), (2,3), (3,4), (4,2), (2,5), (5,0), (6,3), (5,6)])
#g = Graph([(0,1),(0,3), (1,3), (3,3)])
#print(randomWalkBetweenness(g))
#print(criticality(g, 0))

#Ejemplo matriz paper 24
#g = Graph([(0,1), (0,2), (0,3), (1,0), (1,2), (1,3), (2,0), (2,1), (2,3), (2,4), (3,0), (3,1), (3,2), (4,6),(5,4), (6,5), (6,7), (7,1)], directed = True)
#print(entropyRank(g))

#g = Graph([(0,1), (1,2) , (1,3), (3,4)], directed =True)
#g = Graph([(0,1), (1,2), (1,3)])
#m = np.array(g.get_adjacency().data)
#print(mcv(m, [], [],g.is_directed()))




