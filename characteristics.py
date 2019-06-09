from igraph import *
import random
import numpy as np
import math
from scipy import linalg
import itertools
import sys

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

def randomRobustnessIndex(g, m):
    """
    g: Graph
    m:
    return
    """
    v = g.vcount()
    sum = 0
    for i in range(m): #Each random attacks
        for j in range(1, v+1): #Remove j random vertices
            aux = g.copy()
            for k in range(j): #Choose a random vertex j times
                numberV = aux.vcount()
                vertex = random.randint(0, numberV - 1) #Choose a random vertex
                aux.delete_vertices([vertex])
            clusters = aux.components()
            maxCluster = 0
            for c in clusters:
                l = len(c)
                if l > maxCluster:
                    maxCluster = l
            sum += maxCluster / v
    return (1/m) * (1/v) * sum

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
    
    newComponents = aux.components()
    nComponents = []

    #Change vertices numbering to original
    for i in range(len(newComponents)):
        l = []
        for j in range(len(newComponents[i])):
            index =  np.where(num == newComponents[i][j])[0][0]
            l.append(index)
        nComponents.append(l)
    new = 0
    for comp in nComponents:
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
        size = len(c)
        if(size > before):
            before = size
    aux = p(g)
    AComp = aux.components()
    after = 0
    for c in AComp:
        size = len(c)
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
    N0 = 0
    vertices = g.vcount()
    for i in range(vertices): #count number of ordered pairs of vertices
        for j in range(vertices):
            if(i != j):
                if(g.vertex_disjoint_paths(i, j, neighbors = "ignore") != 0): #If there is a path between vertex i and j
                    N0 += 1
    
    aux = g.copy()
    aux.delete_vertices([v])
    Nv = 0
    nvertices = aux.vcount()
    for i in range(nvertices): #count number of ordered pairs of vertices
        for j in range(nvertices):
            if(i != j):
                if(aux.vertex_disjoint_paths(i, j, neighbors = "ignore") != 0): #If there is a path between vertex i and j
                    Nv += 1

    return (N0 - Nv)/ N0

def fragmentation(g, strategy, args):
    """
    g: Graph
    strategy: Function that makes the removal, returns a graph. strategy(g, args)
    args: Arguments for strategy
    return: 
    """
    N = g.vcount()
    if(N == 1):
        return -1 #error
    removed = strategy(g, args)
    clusters = removed.components()
    sum = 0
    for comp in clusters:
        sum+= len(comp)
    return sum / (N * (N - 1))

def selfSufficiency(g, l):
    """
    g: Graph
    l: the set of services available locally  the set of nonlocal services for each vertex. List = [[[A(v_0)], [N(v_0)]], ... , [[A(v_n-1)], [N(v_n-1)]]]
    return:
    """
    comp = g.components()
    for c in comp:
        for vertex in c:
            nonLocal = l[vertex][1] #List of nonlocal services needed at v
            for n in nonLocal:
                cond = False
                for rest in c:
                    if (rest != vertex):
                        if(n in l[rest][0]):
                            cond= True
                            break
                if (not cond):
                    return False
    return True

def kVertexFailureResilience(g, l, k):
    """
    g: Graph
    l: the set of services available locally  the set of nonlocal services for each vertex. List = [[[A(v_0)], [N(v_0)]], ... , [[A(v_n-1)], [N(v_n-1)]]]
    k: Number of vertices that fail
    return:
    """
    if k == 0:
        return selfSufficiency(g, l)
    v = g.vcount()
    if k > v:
        raise Exception('Number of vertices to fail can not be greater than the total vertices')
    nodes = list(range(v))
    for i in range(1, k+1):
        combinations = list(itertools.permutations(nodes, i))
        for comb in combinations:
            auxGraph = g.copy()
            auxList = l.copy()
            auxGraph.delete_vertices(comb)
            for vertex in sorted(comb, reverse= True):
                del auxList[vertex]
            if (not selfSufficiency(auxGraph, auxList)):
                return False
    return True


def resilience(g, l, function):
    """
    g: Graph
    l: the set of services available locally  the set of nonlocal services for each vertex. List = [[[A(v_0)], [N(v_0)]], ... , [[A(v_n-1)], [N(v_n-1)]]]
    function:
    return: The largest k for which g is function resilient
    """
    k = 1
    while(True):
        if(function(g, l, k)):
            k += 1
        else:
            return k - 1

def vertexResilience(g, l):
    """
    g: Graph
    l: the set of services available locally  the set of nonlocal services for each vertex. List = [[[A(v_0)], [N(v_0)]], ... , [[A(v_n-1)], [N(v_n-1)]]]
    return: The largest k for which g is k vertex-failure resilient
    """
    return resilience(g, l, kVertexFailureResilience)

def kEdgeFailureResilience(g, l, k):
    """
    g: Graph
    l: the set of services available locally  the set of nonlocal services for each vertex. List = [[[A(v_0)], [N(v_0)]], ... , [[A(v_n-1)], [N(v_n-1)]]]
    k: Number of edges that fail
    """
    if k == 0:
        return selfSufficiency(g, l)
    e = g.ecount()
    if k > e:
        raise Exception('Number of edges to fail can not be greater than the total edges')
    edges = list(range(e))
    for i in range(1, k+1):
        combinations = list(itertools.permutations(edges, i))
        for comb in combinations:
            auxGraph = g.copy()
            auxList = l.copy()
            for edge in sorted(comb, reverse= True):
                auxGraph.delete_edges(edge)
            if (not selfSufficiency(auxGraph, auxList)):
                return False
    return True

def edgeResilience(g, l):
    """
    g: Graph
    l: the set of services available locally  the set of nonlocal services for each vertex. List = [[[A(v_0)], [N(v_0)]], ... , [[A(v_n-1)], [N(v_n-1)]]]
    return: The largest k for which g is k edge-failure resilient
    """
    return resilience(g, l, kEdgeFailureResilience)

def getSimplePath(g, s, d, seed):
    """
    Auxiliary function
    g: Graph
    s: Source vertex
    d: Destination vertex
    seed: Seed for randomize, if seed = 0, this parameter is not used
    return: [[List of nodes of the path], [List of edges of the path]]
    """
    if(g.vertex_disjoint_paths(s,d, neighbors = "ignore") == 0):
        raise Exception('There is no path between s and d')
    nodes = [s]
    edges = []
    actual = s
    if(seed):
        random.seed(seed)
    while(actual != d):
        neighbors = g.neighbors(actual, mode= "out")
        if (len(neighbors) == 0):
            nodes = [s]
            edges = []
            actual = s
        else:
            new = random.choice(neighbors)
            if not new in nodes:
                edgeId = g.get_eid(actual, new)
                nodes.append(new)
                edges.append(edgeId)
                actual = new
            else:
                nodes = [s]
                edges = []
                actual = s
    return [nodes, edges]

def pathDiversity(g, s, d, seed):
    """
    g: Graph
    s: Source vertex
    d: Destination vertex
    seed: Seed for randomize
    """
    if(g.vertex_disjoint_paths(s,d, neighbors = "ignore") == 0):
        raise Exception('There is no path between s and d')
    L0 = g.get_shortest_paths(s, d, output="epath")[0]
    N0 = g.get_shortest_paths(s, d)[0]
    Nk, Lk = getSimplePath(g, s, d, seed)
    L = list(set(L0) & set(Lk)) #Intersection
    N = list(set(N0) & set(Nk)) #Intersection
    return 1 - (len(L) + len(N))/(len(L0) + len(N0))

def percolatedPath(g, s, d, state):
    """
    g: Graph
    s: Source vertex
    d: Destination vertex
    state: Name of the parameter that indicates the percolation state
    return: The shortest path from s to d such that s is infected
    """
    if(g.vs[s].attributes()[state] == 0):
        raise Exception('Vertex s is not percolated')
    if(g.vertex_disjoint_paths(s, d, neighbors = "ignore") == 0):
        raise Exception('There is no path between s and d')
    return g.get_shortest_paths(s, d)[0]

def percolationCentrality(g, v, state):
    """
    g: Graph
    v: Vertex
    state: Name of the parameter that indicates the percolation state
    """
    n = g.vcount()
    if (n == 2):
        raise Exception('Number of vertices can not be 2')
    xv = g.vs[v].attributes()[state]

    #Sum of x_i
    sumxi = 0
    for i in range(n):
        sumxi += g.vs[i].attributes()[state]
    
    sum = 0
    for s in range(n):
        for r in range(n):
            if(s != v and v!= r):
                paths = g.get_all_shortest_paths(s, r)
                omegasr_v = 0
                for p in paths:
                    if v in p:
                        omegasr_v += 1
                sum += (omegasr_v / len(paths)) * (g.vs[s].attributes()[state] / (sumxi - xv))

    return (1 / (n - 2)) * sum

def getDegreeDistribution(g):
    h = g.degree_distribution() #Degree distribution
    bins = list(h.bins())
    n = g.vcount()
    p_k = np.zeros(n)
    for b in bins:
        min = b[0]
        p_k[math.floor(min)] += b[2]
    
    acump_k = np.add.accumulate(p_k)
    return acump_k / n

def degreeEntropy(g):
    """
    g: Graph
    return:
    """
    pk = getDegreeDistribution(g)
    sum = 0
    for i in range(1, g.vcount()):
        sum += pk[i] * math.log(pk[i])
    return -sum

def relativeEntropy(g):
    n = g.vcount()
    pk = getDegreeDistribution(g)
    sum = 0
    for i in range(n):
        sum += pk[i] * math.log(n * pk[i])
    return sum

def geographicalDiversity(g, p):
    """
    g: Graph
    p: Path between two nodes : p = [L, N], where L =  links id, N = nodes id
    return: The minimun distance between any node members of vector p and the shortest path
    """
    N = p[1]
    s = N[0]
    d = N[len(N) -1]
    shortestPath = g.get_shortest_paths(s, d)[0]
    m = float('inf')
    cond = False
    for node in N:
        for node_s in shortestPath:
            if(g.vertex_disjoint_paths(node, node_s, neighbors = "ignore") != 0): #If exist a path between node and node_s
                distance = g.shortest_paths_dijkstra(node, node_s)[0]
                if(distance < m):
                    m = distance
                    cond = True
    
    if(cond):
        return m
    else:
        raise Exception('There is no path between any node members of vector p and the shortest path')

def getAllSimplePaths(g, s, d):
    """
    g: Graph
    s: Source vertex
    d: Destination vertex
    """
    #TODO
    return []

def effectiveGeographicalPathDiversity(g, s, d, l):
    """
    g: Graph
    s: Source node
    d: Destination node
    l: lambda value
    return:
    """
    k = 0
    simplePaths = getAllSimplePaths(g, s, d)
    for path in simplePaths:
        k += geographicalDiversity(g, path)
    
    return 1 - math.exp(- l * k)

def totalGraphGeographicalDiversity(g, l):
    """
    g: Graph
    l : lambda value for effective geographical path diversity
    return: The average of effective geographical path diversity of all vertex pairs within g.
    """
    nodes = list(range(g.vcount()))
    pairs = list(itertools.permutations(nodes, 2)) #This list does not contain pairs (v,v)
    sum = 0
    for p in pairs:
        sum += effectiveGeographicalPathDiversity(g, p[0], p[1], l)
    return sum / len(pairs)

def compensatedTotalGeographicalGraphDiversity(g, l):
    """
    g: Graph
    l : lambda value for total graph geographical path diversity
    return: The total graph geographical diversity weighted by the total number of edges of g
    """
    e = g.ecount()
    return totalGraphGeographicalDiversity(g, l) * e

def functionality(g, i, v, attack):
    """
    g: Graph
    i: Attack numer i that eliminates vertex v_i
    v: Vertex
    attack: Attack sequence
    """
    v = g.vcount()
    if(i > v):
        raise Exception('i can not be greater than the total number of vertices')
    result = np.zeros(i + 1)
    result[0] = 1
    k = 1
    aux = g.copy()
    aux.delete_vertices([attack[0]])
    while(k <= i):
        distance = aux.shortest_paths_dijkstra(v, k)[0]
        degree = aux.degree(k)
        result[k] = result[k - 1] - (1 / ((distance ** 2) * degree)) * result[k - 1]
        k += 1
        aux.delete_edges([attack[k]])
    return result[i]

def functionalityLoss(g, v, attack):
    """
    g: Graph
    v: Vertex
    k: Number of attacks
    attack: Attack sequence
    return: The sum of the differences of vertex functionality before and after the attack over the sequence of k attacks
    """
    sum = 0
    k = len(attack)
    for i in range(1, k):
        sum += functionality(g, i -1, v, attack) - functionality(g, i , v, attack)
    return sum

def globalFunctionalityLoss(g, attack):
    """
    g: Graph
    attack: Attack sequence
    edge: Condition that indicates if the attack is for edges
    """
    sum = 0
    v = g.vcount()
    for i in range(v):
        if i not in attack:
            sum += functionalityLoss(g, i, attack)
        
def shortestTemporalDistance(g, s, d, t2):
    """
    g: Graph
    s: Source vertex
    d: Destination vertex
    t2: 
    return:
    """
    if(g.vertex_disjoint_paths(s, d, neighbors = "ignore") != 0): #If there is a way between vertex s and d
        path_len = g.shortest_paths_dijkstra(s, d)[0]
        if(path_len > t2):
            return float('inf')
        else:
            return path_len

def temporalEfficiency(g, t2):
    """
    g: Graph
    t2:
    return:
    """
    sum = 0
    v = g.vcount()
    for i in range(v):
        for j in range(v):
            if(i != j):
                if(g.vertex_disjoint_paths(i, j, neighbors = "ignore") != 0):
                    sum += (1 / shortestTemporalDistance(g, i, j, t2))
    return (1 / (v * (v - 1))) * sum

def deltaEfficiency(g, i):
    """
    g: Graph
    i: Vertex
    return: Efficiency of vertex i
    """
    sum = 0
    v = g.vcount()
    for j in range(v):
        if (j != i):
            if(g.vertex_disjoint_paths(i, j, neighbors = "ignore") != 0):
                    sum += (1 / shortestTemporalDistance(g, i, j, float('inf')))
    return (1 / (v - 1)) * sum

def fragility(g, t2):
    """
    g: Graph
    t2:
    return
    """
    sum = 0
    v = g.vcount()
    for k in range(v):
        g_k = g.copy()
        v_k = list(range(v))
        v_k.remove(k)
        #TODO: remove all links connected to node k
        for i in v_k:
            for j in v_k:
                if(i != j and g.vertex_disjoint_paths(i, j, neighbors = "ignore") != 0 and g_k.vertex_disjoint_paths(i, j, neighbors = "ignore") != 0):
                    sum += (1 / shortestTemporalDistance(g, i, j, t2)) - (1 / shortestTemporalDistance(g_k, i, j, t2))
    return (1 / (v * (v - 1) * (v - 2)) * sum) 

def dynamicFragility(g, t2):
    """
    g: Graph
    t2:
    """
    return 1 - fragility(g, t2)

def criticalityOfVertex(g, v, w):
    """
    g: Graph
    v: Vertex
    w: String, name of the atribute for link weight
    return: Criticality of vertex v
    """
    sum = 0
    neighbors = g.neighbors(v, mode= 'IN') #TODO: revisar que sean los incidentes
    for n in neighbors:
        link_id = g.get_eid(n, v)
        sum += g.es[link_id].attributes()[w]
    return g.betweenness(v, directed=g.is_directed(), weights=w) / sum

def criticalityOfEdge(g,i, j, w):
    """
    g: Graph
    i: Source vertex of the edge
    j: Incident vertex
    w: String, name of the atribute for link weight
    return: Criticality of edge e
    """
    link_id = g.get_eid(i,j)
    return edgeBetweenness(g, i, j) / g.es[link_id].attributes()[w]

def edgeBetweenness(g, s, d):
    """
    g: Graph
    s: Source vertex
    d: Destination vertex
    return:
    """
    #TODO
    return 0

def networkCriticality(g, w, onlyEdges=False, both=False):
    """
    g: Graph
    w: String, name of the atribute for link weight
    onlyEdges: If false, then only vertices
    both
    return: The sum of criticalities over all elements of G (vertices, edges or both)
    """
    sum = 0
    v = g.vcount()
    for i in range(v):
        if (not onlyEdges):
            sum += criticalityOfVertex(g, i, w)
        if both:
            for j in range(v):
                if (i != j):
                    sum += criticalityOfEdge(g, i, j, w)
    return sum

    


        














    

            







    

    

       
#g = Graph([(0,1), (0,2), (2,3), (3,4), (4,2), (2,5), (5,0), (6,3), (5,6)])
#g = Graph([(0,1),(0,3), (1,3), (3,3)])
#print(randomWalkBetweenness(g))
#print(criticality(g, 0))

#Ejemplo matriz paper 24
#g = Graph([(0,1), (0,2), (0,3), (1,0), (1,2), (1,3), (2,0), (2,1), (2,3), (2,4), (3,0), (3,1), (3,2), (4,6),(5,4), (6,5), (6,7), (7,1)], directed = True)
#print(entropyRank(g))

#g = Graph([(0,1), (2,1), (0,4), (2,5), (0,3),(5,3),(5,4)], directed=True)
#g = Graph([(0,2), (0,4), (1,5), (2,1), (3,1), (5,3)], directed=True)

g = Graph([(0,1), (0,2), (1,3), (1,2), (2,3), (0,5)])
nodes = list(range(5))
pairs = list(itertools.permutations(nodes, 2)) 
print(pairs[0][0])
print(pairs[0][1])