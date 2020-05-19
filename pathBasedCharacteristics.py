import itertools
import random
from igraph import *
from auxiliaryFunctions import *

def localConnectivity(graph):
    giant_component = getGiantComponent(graph)
    try:
        connectivity_graph = localConnectivityAux(graph)
        connectivity_subgraph = localConnectivityAux(giant_component)
        result = connectivity_graph / connectivity_subgraph
    except ZeroDivisionError:
        result = None
    return result

def globalConnectivity(graph):
    n = graph.vcount()
    try:
        comp = graph.components()
        result = max(comp.sizes()) / n
    except ZeroDivisionError:
        result = None
    return result

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
        return None
    removed = strategy(g, args)
    clusters = removed.components()
    sum = 0
    for comp in clusters:
        sum+= len(comp)
    return sum / (N * (N - 1))

def selfSufficiency(g, l=None): # auxiliary measure
    """
    g: Graph
    l: the set of services available locally  the set of nonlocal services for each vertex. List = [[[A(v_0)], [N(v_0)]], ... , [[A(v_n-1)], [N(v_n-1)]]]
    """
    if l == None:
        l = makeRandomServices(g)
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
        return None
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
        prov_sj = providesSj(l, sj) #List of nodes that provides service dj
        sigma_j = []
        for v in dj:
            alpha_vj = minimumWeightstNodeCutset(v, prov_sj, aux_g)
            sigma_j.append(alpha_vj)
        if len(sigma_j) != 0:
            T.append(np.min(sigma_j))
    if len(T) == 0:
        return None
    return np.min(T) - 1


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
        return None
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
        prov_sj = providesSj(l, sj) #List of nodes that provides service dj
        gamma_j = []
        for v in dj: #v is a node            
            alpha_vj = minimumWeightstEdgeCutset(v, prov_sj, aux_g)
            gamma_j.append(alpha_vj)
        if len(gamma_j) != 0:
            O.append(np.min(gamma_j))
    if len(O) == 0:
        return None
    return np.min(O) - 1

def pathDiversity(g, d=None, s=0,  seed=1):
    """
    g: Graph  
    d: Destination vertex
    s: Source vertex, default=0
    seed: Seed for randomize, default = 0
    """
    if d is None:
        n= g.vcount()
        d = random.randint(0, n-1)

    if(g.vertex_disjoint_paths(s,d, neighbors = "ignore") == 0):
        return None
          
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

def percentageOfNoncriticalNodes(CN, g2, interdep= None, p=0.85, admittance = None, voltage = None, tolerance = 0.5, a_t = 0.2, a_p = 0.2):
    """
    Note: CN and PG must be the same number of nodes
    CN: Communication network
    g2: Power grid network
    interdep: List of interdependencies, format [(CN, PG)], each node can have at most one interdependent link
    p: Percentage of interdependent nodes, only used when interdep is None
    admittance: Edge attribute for admittance
    voltage: Vertex attribute for voltage
    generation_nodes: List of generation nodes
    tolerance: Tolerance threshold
    a_t: Tolerance parameter of communication network 
    a_p: Tolerance parameter of nodes in power grid network

    """
    PG = g2
    #Check graphs
    if CN.vcount() != PG.vcount():
        return None

    # Set all variables for PG
    if admittance is None:
        admittance = 'admittance'
        PG = generateWeight(PG, edge = True, vertex = False, name = admittance)   
    
    if voltage is None:
        numberOfGenerationNodes = int(PG.vcount() / 10) #10% of the nodes
        nodes = list(range(PG.vcount()))
        generation_nodes = random.sample(nodes, numberOfGenerationNodes)
        # Set positive voltage for generation nodes
        voltage = 'voltage'
        PG = generateVoltages(PG, voltage, generation_nodes)

    # Set id for each node in both networks
    PG = setId(PG, name='id')
    CN = setId(CN, name='id')

    if interdep is None:
        CN_dict = {}
        PG_dict = {}
        # Set interdependencies
        n_interped = int(p * CN.vcount()) # Number of interdependent nodes
        CN_nodes = list(range(CN.vcount()))
        PG_nodes = list(range(PG.vcount()))

        CN_interp_nodes = random.sample(CN_nodes, n_interped)
        PG_interp_nodes = random.sample(PG_nodes, n_interped)
        for CN_node in CN_interp_nodes:
            PG_node = random.choice(PG_interp_nodes)
            PG_interp_nodes.remove(PG_node)
            CN_dict[CN_node] = PG_node
            PG_dict[PG_node] = CN_node
    else:
        CN_dict, PG_dict = toDict(interdep)

    # Maximum capacities of each network
    CN_initial_load = communicationNetworkLoad(CN)
    CN_max_capacity = (1 + a_t) * CN_initial_load
    
    PG_initial_load = powerGridNetworkLoad(PG, voltage, admittance)
    PG_max_capacity = (1 + a_p) * PG_initial_load

    CN_n = CN.vcount()
    PG_n = PG.vcount()
    n = PG_n + CN_n

    acc = 0
    for v in range(CN_n):
        # Delete node and counterpart in PG
        PG_aux = PG.copy()
        CN_aux = CN.copy()
        n_v = 0
        if v in CN_dict:
            n_v = 1
            PG_neighbor = CN_dict[v]
            PG_aux.delete_vertices(PG_neighbor)
        CN_aux.delete_vertices(v)

        n_v += cascadingFailures(CN_aux, PG_aux, CN_max_capacity, PG_max_capacity, voltage, admittance, CN_dict, PG_dict)
        if n_v / n < tolerance:
            acc += 1
    
    for v in range(PG_n):
        # Delete node and counterpart in CN
        PG_aux = PG.copy()
        CN_aux = CN.copy()
        n_v = 0
        if v in PG_dict:
            n_v = 1
            CN_neighbor = PG_dict[v]
            CN_aux.delete_vertices(CN_neighbor)
        PG_aux.delete_vertices(v)

        n_v += cascadingFailures(CN_aux, PG_aux, CN_max_capacity, PG_max_capacity, voltage, admittance, CN_dict, PG_dict)
        if n_v / n < tolerance:
            acc += 1
    res = acc / n
    return res

def treeness(g):
    """
    g: Graph
    """
    Gc = getDAG(g)
    Gc_r = getDAG(g, reverse=True)
    hf = H_f(Gc)
    hb = H_f(Gc_r)
    m = max(hf, hb)
    if m == float('inf'):
        return 0
    try:
        return (hf - hb) / m
    except ZeroDivisionError:
        return None
