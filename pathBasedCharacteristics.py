import itertools
import random
from igraph import *
from auxiliaryFunctions import *

def pairwiseDisconnectivityIndex(g, v=0):
    """
    g: Graph
    v: Vertex, default = 0
    return: The pairwise disconnectivity of vertex v
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

def fragmentation(g, strategy=perturbationFunction, args=1):
    """
    g: Graph
    strategy: Function that makes the removal, returns a graph. strategy(g, args)
    args: Arguments for strategy
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

def selfSufficiency(g, l=0):
    """
    g: Graph
    l: the set of services available locally  the set of nonlocal services for each vertex. List = [[[A(v_0)], [N(v_0)]], ... , [[A(v_n-1)], [N(v_n-1)]]]
    """
    if l == 0:
        l = makeEmptyServices(g)
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

def kVertexFailureResilience(g, l=0, k= 1):
    """
    g: Graph
    l: the set of services available locally  the set of nonlocal services for each vertex. List = [[[A(v_0)], [N(v_0)]], ... , [[A(v_n-1)], [N(v_n-1)]]] 
    k: Number of vertices that fail, default = 1
    """
    if l == 0:
        l = makeEmptyServices(g)
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

def vertexResilience(g, l=0):
    """
    g: Graph
    l: the set of services available locally  the set of nonlocal services for each vertex. List = [[[A(v_0)], [N(v_0)]], ... , [[A(v_n-1)], [N(v_n-1)]]]
    return: The largest k for which g is k vertex-failure resilient
    """
    if l == 0:
        l = makeEmptyServices(g)
    return resilience(g, l, kVertexFailureResilience)

def kEdgeFailureResilience(g, l=0, k=1):
    """
    g: Graph
    l: the set of services available locally  the set of nonlocal services for each vertex. List = [[[A(v_0)], [N(v_0)]], ... , [[A(v_n-1)], [N(v_n-1)]]]
    k: Number of edges that fail, default = 1
    """
    if l == 0:
        l = makeEmptyServices(g)
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

def edgeResilience(g, l=0):
    """
    g: Graph
    l: the set of services available locally  the set of nonlocal services for each vertex. List = [[[A(v_0)], [N(v_0)]], ... , [[A(v_n-1)], [N(v_n-1)]]]
    return: The largest k for which g is k edge-failure resilient
    """
    if l == 0:
        l = makeEmptyServices(g)
    return resilience(g, l, kEdgeFailureResilience)

def pathDiversity(g, d= 1, s=0,  seed=1):
    """
    g: Graph  
    d: Destination vertex
    s: Source vertex, default=0
    seed: Seed for randomize, default = 0
    """
    if(g.vertex_disjoint_paths(s,d, neighbors = "ignore") == 0):
        raise Exception('There is no path between s and d')
    L0 = g.get_shortest_paths(s, d, output="epath")[0]
    N0 = g.get_shortest_paths(s, d)[0]
    Nk, Lk = getSimplePath(g, s, d, seed)
    L = list(set(L0) & set(Lk)) #Intersection
    N = list(set(N0) & set(Nk)) #Intersection
    return 1 - (len(L) + len(N))/(len(L0) + len(N0))

def percolatedPath(g, d=1, s=0, state='state'):
    """
    g: Graph
    d: Destination vertex
    s: Source vertex, default = 0
    state: Name of the parameter that indicates the percolation state, default='state'
    return: The shortest path from s to d such that s is infected
    """
    if(g.vs[s].attributes()[state] == 0):
        raise Exception('Vertex s is not percolated')
    if(g.vertex_disjoint_paths(s, d, neighbors = "ignore") == 0):
        raise Exception('There is no path between s and d')
    return g.get_shortest_paths(s, d)[0]

def percolationCentrality(g, v=0, state='state'):
    """
    g: Graph
    v: Vertex, default=0
    state: Name of the parameter that indicates the percolation state, default='state'
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