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

def splittingNumber(g, k, dev= 1):
    """
    g: Graph
    k:
    dev:
    return: The average number of edges that need to be removed to break  g into k connected components
    """

    numOfLoops = 0 #TODO: definir
    givenComponents = k

    i = 0
    numOfComponents = 0
    cond = True
    l = g.ecount()
    k = 3 #TODO: Calculate
    numOfLinks = np.zeros(k)
    while (cond):
        auxGraph = g.copy() #make all links functioning
        numOfComponents = len(auxGraph.components()) #Calculate numOfcomponents
        numOfLinks[i] = 0

        while (numOfComponents < givenComponents):
            #Choose randomly an uninterrupted link
            numberTotalLinks = auxGraph.ecount()
            linkId = random.randint(0, numberTotalLinks - 1)            

            numOfLinks[i] += 1
            auxGraph.delete_edges([linkId]) #interrup the link
            numOfComponents = len(auxGraph.components())#compute numOfComponents

            if(i > 1 and i > numOfLoops):
                mean1 = np.mean(numOfLinks[1:i+1]) #mean value of numOfLinnks[1:i]
                mean2 = np.mean(numOfLinks[1:i]) #mean value of numOfLinnks[1:i - 1]
                if(abs(mean1 - mean2) < dev):
                    return (mean1 + mean2)/2
        i += 1

def robustnessMeasure53(g):
    """
    g: Graph
    return
    """
    aux = g.copy() #The function doesn't modify the graph
    n = aux.vcount()
    if (n < 1):
        return -1 #Error
    Ct = 0
    for i in range(n):
        degrees = aux.degree()
        v = np.argmax(degrees) #Vertex of maximun degree
        aux.delete_vertices([v])
        C = aux.components() #Components of the graph
        Ci = 0 #Order of the biggest component
        for comp in C:
            compOrder = len(comp)
            if(compOrder > Ci):
                Ci = compOrder 
        Ct += Ci
    return Ct / n

def connectivityRobustnessFunction(g, k):
    """
    g: Graph
    k: Number of vertices removed
    return 
    """
    n = g.vcount()
    if(k > n):
        return -1 #Error
    aux = g.copy() #The function doesn't modify the original graph
    for i in range(k): #Delete k random vertices
        numberV = aux.vcount()
        v = random.randint(0, numberV - 1) #Choose a random vertex
        aux.delete_vertices([v]) #Delete vertex
    C = aux.components() #Components of the graph
    S = 0
    for comp in C: #Get the largest connected component
        compOrder = len(comp)
        if(compOrder > S):
                S = compOrder
    return S / (n - k)

def kResilienceFactor(g, k):
    """
    g: Graph
    k:
    return: The percentage of connected components of g that remain connected after the removal k - 1 vertices 
    """
    n = g.vcount()
    if(k - 1 > n):
        return -1 #Error
    aux = g.copy() #The function doesn't modify the original graph
    originalComponents = g.components()
    original = len(originalComponents) #Number of original connected components
    num = np.arange(n) #New numbering after removing vertices
    for i in range(k - 1): #Remove k - 1 vertices
        numberV = aux.vcount()
        v = random.randint(0, numberV -1)
        aux.delete_vertices([v]) #Delete vertex
        #Change vertices numbering
        index = np.where(num == v)[0][0]
        num[index] = -1
        for j in range(index, len(num)):
            if(num[j] != -1):
                num[j] -= 1
    
    newComponents = aux.components() #TODO: Modificar para tener vertices con números originales
    #Change vertices numbering to original
    for i in range(len(newComponents)):
        for j in range(len(newComponents[i])):
            index =  np.where(num == newComponents[i][j])[0][0]
            newComponents[i][j] = np.where(num == newComponents[i][j])
    new = 0
    for comp in newComponents:
        if comp in originalComponents:
            new += 1
    return (new/original) * 100

def resilienceFactor(g):
    n = g.vcount()
    if( n < 3):
        return -1 #Error
    result = np.zeros(n-2)
    for i in range(2, n):
        auxGraph = g.copy()
        result[i-2] = kResilienceFactor(auxGraph, i)
    return np.mean(result)

def pertubationScore(g, p):
    """
    g: Graph
    p: Perturbation function
    return:
    """
    BComp = g.components()
    before = 0
    for c in BComp:
        size = c.vcount()
        if(size > before):
            before = size
    aux = p(g)
    AComp = aux.components()
    after = 0
    for c in AComp:
        size = c.vcount()
        if(size > after):
            after = size
    return (before - after)/ before

def preferentialPerturbation(g1, g2):
    """
    g1: Graph
    g2: Graph
    return:
    """
    return

def maximunPerturbationScore(g1,g2):
    """
    g1: Graph
    g2: Graph
    return:
    """
    return

def pairwiseDisconnectivityIndex(g, v):
    """
    g: Graph
    v: Vertex
    return
    """
    

    

       
#g = Graph([(0,1), (0,2), (2,3), (3,4), (4,2), (2,5), (5,0), (6,3), (5,6)])
#g = Graph([(0,1),(0,3), (1,3), (3,3)])
#print(randomWalkBetweenness(g))
#print(criticality(g, 0))

#Ejemplo matriz paper 24
#g = Graph([(0,1), (0,2), (0,3), (1,0), (1,2), (1,3), (2,0), (2,1), (2,3), (2,4), (3,0), (3,1), (3,2), (4,6),(5,4), (6,5), (6,7), (7,1)], directed = True)
#print(entropyRank(g))

#g = Graph([(0,1), (2,1), (0,4), (2,5), (0,3),(5,3),(5,4)], directed=True)
g = Graph([(1,3)])
print(kResilienceFactor(g, 3))
