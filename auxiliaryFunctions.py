import random
import numpy as np
from scipy import linalg
from igraph import *
from itertools import combinations

def randomWalk(g, i=0, t=1, s= 0):
    """
    g: Graph
    i: Source vertex, default=0
    t: Target, default=1
    s: Seed for randomize, default = 0
    return: Random walk between s and t, list where each element is a vertex
    """
    if (g.vertex_disjoint_paths(i, t, neighbors = "ignore") == 0): #If there is no path between s and t
        return []
    if(s):
        random.seed(s)
    l = [i]
    actual = i
    while(actual != t):
        neighbors = g.neighbors(actual, mode='OUT')
        if len(neighbors) == 0:
            actual = i
            l = [i]
        else:
            n = random.choice(neighbors)
            l.append(n)
            actual = n
    return l

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

def allCovers(g):
    r = []
    p = []
    covers(g, p, r, g.is_directed())
    return r

def covers(g, p=[], r=[], directed=False):
    """
    Get all vertex covers of graph g, the covers are stored in r
    g: Graph
    p : Partial cover
    r: Variable where the covers are stored
    directed: Boolean that indicates if the graph is directed
    """
    e = g.ecount()
    if (e == 0):
        p.sort()
        if(p not in r):
            r.append(p)
    else:
        for edge in g.es:
            u = edge.tuple[0]
            v = edge.tuple[1]
            g_aux = g.copy()
            g_aux.delete_edges([(u,v)])
            if not directed:
                if (u not in p):
                    n_u = g_aux.neighbors(u)
                    g_u = g_aux.copy()
                    for n in n_u:
                        g_u.delete_edges([(u, n)])
                    p_u = p.copy()
                    p_u.append(u)
                    covers(g_u, p_u, r, directed)
            if (v not in p):
                n_v = g_aux.neighbors(v, mode="IN")
                g_v = g_aux.copy()
                for n in n_v:
                    g_v.delete_edges([(n, v)])
                p_v = p.copy()
                p_v.append(v)
                covers(g_v, p_v, r, directed)

def mcv(g):
    """
    g: Graph
    return: The set of minimal vertex covers of G
    """
    covers = allCovers(g)
    r = []
    for c_cover in covers:
        covers_aux = covers.copy()
        #Delete current cover from cover auxiliary list
        covers_aux.remove(c_cover)
        cond = False #Does c_cover contain a cover
        for cover in covers_aux:
            cov_cond = all(elem in c_cover for elem in cover) #Is the cover contain in c_cover
            if cov_cond:
                cond = True
        if(not cond):
            r.append(c_cover)
    if not r:
        return None
    return r

def MCV(g):
    """
    g: Graph
    return: The set of minimum vertex covers of G
    """
    mcv_covers = mcv(g)
    result = []
    min = len(mcv_covers[0])
    for cover in mcv_covers:
        l = len(cover)
        if l < min:
            min = l
    for cover in mcv_covers:
        if(len(cover) == min):
            result.append(cover)
    if not result:
        return None
    return result

def weightedMatrix(g, w):
    """
    g: Graph
    w: Strings that represent the name of the edge weight attribute in the graph
    return: The weighted adjacency matrix of graph g
    """
    m = g.get_adjacency().data
    n = len(m)
    for i in range(n):
        for  j in range(n):
            if(m[i][j] != 0):
                e_id = g.get_eid(i,j)
                e_weight = g.es[e_id].attributes()[w]
                m[i][j] = e_weight
    return m

def sizeMaxComponent(g):
    """
    g: Graph
    return: The size of the biggest component in g
    """
    comp = g.components()
    size = 0
    for c in comp:
        actual_size = len(c)
        if(actual_size > size):
            size = actual_size
    return size

def perturbationScoreTwo(g1, g2):
    """
    g1: Graph before the perturbation
    g2: Graph after the perturbation
    return: The perturbation score given two graphs
    """
    before = sizeMaxComponent(g1)
    after = sizeMaxComponent(g2)
    return (before - after)/ before

def perturbationFunction(g, seed=1):
    """
    This function perturbates a graph by deleting some vertices
    g: Graph
    return: A perturbated graph
    """
    aux = g.copy()
    v = aux.vcount()
    if(seed):
        random.seed(seed)
    c = random.randint(0,v)
    for i in range(c):
        v_id = random.randint(0, aux.vcount() - 1)
        aux.delete_vertices(v_id)
    return aux

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

def resilience(g, l, function):
    """
    g: Graph
    l: the set of services available locally  the set of nonlocal services for each vertex. List = [[[A(v_0)], [N(v_0)]], ... , [[A(v_n-1)], [N(v_n-1)]]]
    function: kVertexFailureResilience or kEdgeFailureResilience
    return: The largest k for which g is function resilient
    """
    k_i = 1
    while(True):
        print(k_i)
        if(function(g, l, k= k_i)):
            k_i += 1
        else:
            return k_i - 1

def getServices(l):
    """
    List of services of each vertex
    return: List of all services
    """
    s = set([])
    for n in l:
        for a in n[0]:
            s.add(a)
        for na in n[1]:
            s.add(na)
    return list(s)

def computedj(l, j):
    """
    l: List of services of each vertex
    j: Service
    return: List of nodes that need service j
    """
    dj = []
    for i in range(len(l)):
        if j in l[i][1]:
            dj.append(i)
    return dj 

def auxGraphj(g, l, j, node_w = False):
    """
    g: Grpah
    l: List of services
    j : Service
    return: Auxiliary graph_j
    """
    v = g.vcount()
    new_g = Graph()
    s = 0
    #Compute pj
    pj = []
    for n in range(len(l)):
        if(j in l[n][0]):
            pj.append(n)
            

    #Compute Vj
    if (v == 0):
        new_g.add_vertices(1)
    else:
        new_g.add_vertices(v)
        s = v
    #Compute Ej1
    ej1 = g.get_edgelist()
    x_j = list(combinations(pj,2))
    #Delete x_j from edge_list
    for p in x_j:
        ej1.remove(p)
    weights = []
    if not node_w:
        weights = np.ones(len(ej1))
    #Compute Ej2
    ej2 = []
    for n in pj:
        ej2.append((s, n))
        if not node_w:
            weights.append(sys.maxsize)
    new_g.add_edges(ej1)
    new_g.add_edges(ej2)
    if not node_w:
        new_g.es['weight'] = weights
    else:
        weights = np.ones(v)
        weights.append(sys.maxsize)
        new_g.vs['weight'] = weights

    return new_g
    

def getProbabilityDegree(g, k, myMode='ALL'):
    """
    g: Graph
    k: 
    myMode: degree distribution mode, it can be 'ALL', 'IN' or 'OUT'
    return: The probability p(k) of the degree distribution
    """
    h = g.degree_distribution(mode= myMode) #Degree distribution
    bins = list(h.bins())
    acc = 0
    for b in bins:
        min = math.floor(b[0])
        if(min <= k):
            acc += b[2]
        else:
            break   
    return acc / g.vcount()

def getDegreeDistribution(g, myMode='ALL'):
    """
    g: Graph
    myMode: degree distribution mode, it can be 'ALL', 'IN' or 'OUT'
    return: List with the degree distribution
    """
    h = g.degree_distribution(mode= myMode) #Degree distribution
    bins = list(h.bins())
    p_k = np.zeros(len(bins))
    i = 0
    for b in bins:
        p_k[i] += b[2]
        i += 1   
    return p_k / g.vcount()

def maxDegreeProduct(g):
    """
    g: Graph
    return: Edge id of the edge with the max degree product
    """
    m = g.ecount()
    degree = np.zeros(m)
    for i in range(m):
        s = g.get_edgelist()[i: i+1][0][0] #Source vertex of the edge
        d = g.get_edgelist()[i: i+1][0][1] #Destination vertex of the edge
        degree[i] = g.degree(s) * g.degree(d)

    return np.argmax(degree)
    

def maxEdgeBetweenness(g):
    """
    g: Graph
    return: Edge id of the edge with the max edge betweenness
    """
    bet = g.edge_betweenness(directed=g.is_directed())
    return bet.index(max(bet))

def shortestTemporalDistance(g, s, d, t2):
    """
    g: Graph
    s: Source vertex
    d: Destination vertex
    t2: Limit of the time window
    return: the smallest length among all the temporal paths between nodes s and d in time window [0, t2]
    """
    if(g.vertex_disjoint_paths(s, d, neighbors = "ignore") != 0): #If there is a way between vertex s and d
        path_len = g.shortest_paths_dijkstra(s, d)[0][0]
        if(path_len > t2):
            return float('inf')
        else:
            return path_len

def getAllSimplePaths(g, s, d, visited, partial = [], result= []):
    """
    g: Graph
    s: Source vertex
    d: Destination vertex
    visited: List of booleans each represent if the vertex is visited, type: numpy array
    partial: Partial path
    result: actual result
    return: A list of all the simple paths between vertex s and d
    """
    partial.append(s)
    visited[s] = True
    if(s == d):
        result.append(partial)
        return result
    neighbors = g.neighbors(s)
    for n in neighbors:
        partial_aux = partial.copy()
        visited_aux = visited.copy()
        if (not visited[n]):
            result = getAllSimplePaths(g, n, d, visited_aux, partial_aux, result)
            visited[n] = True
    return result

def get_edges(g, k):
    """
    g: Graph
    v: Vertex
    return: A list of all links id connected to vertex k
    """
    neighbors = g.neighbors(k)
    result = []
    for n in neighbors:
        edge_id = g.get_eid(k, n)
        result.append(edge_id)
    return result

def criticalityOfEdge(g, i= 0, j= 1, w='weight'):
    """
    g: Graph
    i: Source vertex of the edge, default = 0
    j: Incident vertex, default = 1
    w: String, name of the atribute for link weight, default = 'weight'
    return: Criticality of edge e, it assumes that edge exists, default edge (0, 1)
    """
    link_id = g.get_eid(i,j)
    edgeBet = g.edge_betweenness(directed=g.is_directed())[link_id]
    return edgeBet / g.es[link_id].attributes()[w]

def criticalityOfVertex(g, v= 0, w= 'weight'):
    """
    g: Graph
    v: Vertex, default = 0
    w: String, name of the atribute for link weight, default = 'weight'
    return: Criticality of vertex v
    """
    sum = 0
    neighbors = g.neighbors(v, mode= 'IN')
    for n in neighbors:
        link_id = g.get_eid(n, v)
        sum += g.es[link_id].attributes()[w]
    return g.betweenness(v, directed=g.is_directed(), weights=w) / sum

def weightFunction(u):
    """
    Auxiliary function for relative area index
    """
    return u ** 2

def maxFlow(v, u):
    """
    Auxiliary function for relative area index
    """
    return v *2 +  u*2

def sortEigenValuesVectors(eigenvalues, eigenvectors, asc=False,absValue= True):
    """
    Auxiliary function
    eigenvalues: Array of eigenvalues
    eigenvectors: Array of eigenvector
    asc: Sort reverse
    abs: Take tha absolute value to compare
    return: Sorted eigen values and eigen vectors, eigenvectors[:,i] is the eigenvector corresponding to the eigenvalues[i]
    """
    if absValue:
        eigenvalues = np.abs(eigenvalues)
    
    eigenvalues = list(eigenvalues)
    eigenvectors = list(eigenvectors)

    pairs = list(zip(eigenvalues, eigenvectors))
    pairs = sorted(pairs, reverse=asc, key=lambda x: x[0])
    unzip = list(zip(*pairs))
    eigenvalues = np.array(unzip[0])
    eigenvectors = np.array(unzip[1])
    return eigenvalues, eigenvectors

def naturalConAux(m):
    """
    Auxiliary function for natural connectivity, used to calculate the natural connectivity of the components of a graph
    m: Adjacency matrix of a component
    """
    eigenvalues = np.linalg.eig(m)[0]
    n = len(eigenvalues)
    exp = np.exp(eigenvalues)
    sum = np.sum(exp)
    return np.log(sum / n)

def changeMatrix(m, r):
    np.sort(r)
    for i in reversed(r):
        np.delete(m, i, 0)
        np.delete(m, i, 1)
    return m
    

def probis(g, i= 0, s= 0, m=[]):
    """
    g: Graph
    i: Vertex, default = 0
    s: State, default = 0
    m: Marked vertices
    return: probability that node i is infected at steady state s.
    """
    neighbors = g.neighbors(i)
    m.append(i)
    acc = 0
    for j in neighbors:
        if j in m:
            acc += getState(g, j, s)
        else:
            acc += probis(g, j, s)
    if(s + acc == 0):
        return 1
    return acc / (s + acc)

def getState(g, v, s):
    """
    g: Graph
    v: Vertex
    s: State
    return: 1 if vertex v is infected at state s, 0 otherwise
    """
    n = g.vcount()
    g_state = ''.join(reversed([str((s >> i) & 1) for i in range(n)]))
    return int(g_state[n - 1 - v])

def y(g, s=0):
    """
    g: Graph
    s: State, default = 0
    return: Fraction of infected nodes at state s
    """
    acc = 0
    v = g.vcount()
    if v == 0:
        return None
    for i in range(v):
        acc += probis(g,i, s)
    return acc / v

def makeEmptyServices(g):
    """
    g: Graph
    return: A list of services per node in a graph [[[], []], ... , [[], []]]
    """
    l = []
    v = g.vcount()
    for i in range(v):
        l.append([[], []])
    return l

def centralityFunction(M):
    """
    Example of centrality function
    M: Adjacency matrix
    """
    n = len(M)
    res = np.zeros(n)
    for i in range(n):
        res[i] = M[i][i] ** 2
    return res