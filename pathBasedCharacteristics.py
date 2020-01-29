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

def selfSufficiency(g, l=None):
    """
    g: Graph
    l: the set of services available locally  the set of nonlocal services for each vertex. List = [[[A(v_0)], [N(v_0)]], ... , [[A(v_n-1)], [N(v_n-1)]]]
    """
    if l == None:
        l = makeRandomServices(g)
        print(l)
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
                    return False, l
    return True, l

def kVertexFailureResilience(g, l=None, k= None):
    """
    g: Graph
    l: the set of services available locally  the set of nonlocal services for each vertex. List = [[[A(v_0)], [N(v_0)]], ... , [[A(v_n-1)], [N(v_n-1)]]] 
    k: Number of vertices that fail, default = random
    """
    if l is None:
        l = makeRandomServices(g)
    if k is None:
        k = random.randint(0, g.vcount())
    if k == 0:
        return selfSufficiency(g, l)
    v = g.vcount()
    if k > v:
        return None, l
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
                return False, l
    return True, l

def vertexResilience(g, l=None):
    """
    g: Graph
    l: the set of services available locally  the set of nonlocal services for each vertex. List = [[[A(v_0)], [N(v_0)]], ... , [[A(v_n-1)], [N(v_n-1)]]]
    return: The largest k for which g is k vertex-failure resilient
    """
    if l == None:
        l = makeRandomServices(g)
    s = getServices(l)
    T = []
    for sj in s:
        aux_g = auxGraphj(g, l, sj , True)
        dj = computedj(l,sj) #List of nodes that needs service sj
        Tj = []
        for v in dj:
            n_j = providesSj(l, sj)
            n_j = onlyReachableNodes(aux_g, v, n_j)
            result = []
            getNodeCutset(v, n_j, aux_g, result, partial=0)
            gamma_vj = np.min(result)
            Tj.append(gamma_vj)
        T.append(np.min(Tj))
    if len(T) == 0:
        return None, l
    return np.min(T) - 1, l


def kEdgeFailureResilience(g, l=None, k=None):
    """
    g: Graph
    l: the set of services available locally  the set of nonlocal services for each vertex. List = [[[A(v_0)], [N(v_0)]], ... , [[A(v_n-1)], [N(v_n-1)]]]
    k: Number of edges that fail, default = random
    """
    if l == None:
        l = makeRandomServices(g)
    if k is None:
        k = random.randint(0, g.ecount())
    if k == 0:
        return selfSufficiency(g, l)
    e = g.ecount()
    if k > e:
        return None, l
    edges = list(range(e))
    for i in range(1, k+1):
        combinations = list(itertools.permutations(edges, i))
        for comb in combinations:
            auxGraph = g.copy()
            auxList = l.copy()
            for edge in sorted(comb, reverse= True):
                auxGraph.delete_edges(edge)
            if (not selfSufficiency(auxGraph, auxList)):
                return False, l
    return True, l

def edgeResilience(g, l=None):
    """
    g: Graph
    l: the set of services available locally  the set of nonlocal services for each vertex. List = [[[A(v_0)], [N(v_0)]], ... , [[A(v_n-1)], [N(v_n-1)]]]
    return: The largest k for which g is k edge-failure resilient
    """
    if l == None:
        l = makeRandomServices(g)
    
    s = getServices(l)
    O = []
    for sj in s:
        aux_g = auxGraphj(g, l, sj)
        dj = computedj(l, sj) #List of nodes that needs service sj
        Oj = []
        for v in dj: #v is a node
            n_j = providesSj(l, sj) #List of nodes that provides service dj
            n_j = onlyReachableNodes(aux_g, v, n_j)
            result = []
            getEdgeCutset(v, n_j, aux_g, result, partial = 0)
            alpha_vj = np.min(result)
            Oj.append(alpha_vj)
        O.append(np.min(Oj))
    if len(O) == 0:
        return None, l
    return np.min(O) - 1, l


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

def percolatedPath(g, s=False, d=False, state=False):
    """
    g: Graph
    s: Source vertex, default = False, if false takes a random vertex
    d: Destination vertex, default = False, if false takes a random vertex
    state: Name of the parameter that indicates the percolation state, default='state'
    return: The shortest path from s to d such that s is infected
    """
    if not state:
        #Set states 0 <= e_i <= 1
        g = setStates(g)
        state = 'state'
    if not s or not d:
        #Choose random s and d
        v = g.vcount()
        s = random.randrange(v)
        vertices = list(range(v))
        vertices.remove(s)
        d = random.choice(vertices)
        #Set s to percolated
        g.vs[s][state] = 1

    if(g.vs[s].attributes()[state] == 0):
        raise Exception('Vertex s is not percolated')
    if(g.vertex_disjoint_paths(s, d, neighbors = "ignore") == 0):
        raise Exception('There is no path between s and d')
    return len(g.get_shortest_paths(s, d)[0])

def percolationCentrality(g, v=False, state=False):
    """
    g: Graph
    v: Vertex, default=False, if false takes a random vertex
    state: Name of the parameter that indicates the percolation state, default=False
    """
    if not state:
        g = setStates(g)
        state= 'state'
    n = g.vcount()
    if not v:
        v = random.randrange(n)
        #TODO: set v to percolated
        g.vs[v][state] = 1
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