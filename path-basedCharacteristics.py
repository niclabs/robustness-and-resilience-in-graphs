import itertools
import random
from igraph import *

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