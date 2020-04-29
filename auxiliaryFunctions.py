import random
import numpy as np
from scipy import linalg
from igraph import *
from itertools import combinations
import networkx as nx

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

def allCovers(g):
    #Complement graph
    complement = Graph.Full(g.vcount())
    edge_list = g.get_edgelist()
    complement.delete_edges(edge_list)

    #Get cliques
    cliques = complement.cliques()

    #Complement of cliques are covers
    covers=[]
    for c in cliques:
        vertex_list = list(range(g.vcount()))
        c = list(c)
        for v in c:
            vertex_list.remove(v)
        covers.append(vertex_list)

    return covers

def MCV(g, mcv_covers):
    """
    g: Graph
    return: The set of minimum vertex covers of G
    """

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
    node_w : Weights at nodes, if False, returns a edge weighted graph
    return: Auxiliary graph_j
    """
    v = g.vcount()
    new_g = Graph()
    s = 0        
    #Compute Vj
    new_g.add_vertices(v +1)
    s = v
    #Compute pj
    pj = providesSj(l, j) 
    #Compute Ej1
    ej1 = g.get_edgelist()
    x_j = list(combinations(pj,2))
    #Delete x_j from edge_list
    for p in x_j:
        if p in ej1:
            ej1.remove(p)
    weights = np.array([])
    if not node_w:
        weights = np.ones(len(ej1))
    #Compute Ej2
    ej2 = []
    for n in pj:
        ej2.append((s, n))
        if not node_w:
            weights = np.append(weights, sys.maxsize)

    new_g.add_edges(ej1)
    new_g.add_edges(ej2) 
    if not node_w:
        new_g.es['weight'] = weights
    else:
        weights = np.ones(v + 1)
        new_g.vs['weight'] = weights

    return new_g

def providesSj(l, j):
    """
    l: List of services of each vertex
    j: Service sj
    return: List of vertices that provides service sj
    """
    result = []
    for v in range(len(l)):
        if j in l[v][0]:
            result.append(v)
    return result

def noPath(g, v, l):
    """
    g: Graph
    v: Vertex
    l: List of nodes
    return: True if there is no path between v and all of the nodes in l
    """
    result = True
    for node in l:
        if(g.vertex_disjoint_paths(v, node, neighbors = "ignore") == 0):
            result = False
    return result

def onlyReachableNodes(g, v, l):
    """
    g: Graph
    v: Vertex
    l: List of nodes
    return: List of nodes in l that are reachable for v
    """
    result = []
    for n in l:
        if(g.vertex_disjoint_paths(v, n, neighbors = "ignore") != 0):
            result.append(n)
    return result

def minimumWeightstNodeCutset(s, nodes, g):
    """
    s: Vertex
    nodes: List of nodes
    g: Graph
    return: Minumin weight node cutset between v and the noeds in 'nodes'
    """
    v = g.vcount()
    all_nodes = list(range(v))
    graph = g.copy()
    #Delete edges between s and nodes in 'nodes'
    for node in nodes:
        if graph.are_connected(s, node):
            graph.delete_edges([(s, node)])
    #Compute minumum weight node cutset
    for i in range(1, v): #For each length of set of vertices
        groups = list(combinations(all_nodes, i)) #Make all posible groups of length i
        for group in groups:
            group = list(group)
            #Valid group
            if s in group:
                continue
            valid_group = True
            for node in nodes:
                valid_group = valid_group and node not in group
            if not valid_group:
                continue
            #Check if group is a cutset
            remain_nodes = all_nodes.copy()
            aux_graph = graph.copy()
            aux_graph.delete_vertices(group)
            for v in group:
                remain_nodes.remove(v)
            s_index = remain_nodes.index(s)
            cut_set = True
            for node in nodes:            
                node_index = remain_nodes.index(node)
                cut_set = cut_set and aux_graph.vertex_disjoint_paths(s_index, node_index, neighbors = "ignore") == 0 #There is no path between s and node               
            if cut_set:
                print(group)
                return i
    return v

def minimumWeightstEdgeCutset(s, nodes, g):
    """
    s: Vertex
    nodes: List of nodes
    g: Graph
    return: Minumin weight node cutset between v and the noeds in 'nodes'
    """
    e = g.ecount()
    all_edges = g.get_edgelist()
    #Compute minumum weight node cutset
    for i in range(1, e): #For each length of set of vertices
        groups = list(combinations(all_edges, i)) #Make all posible groups of length i
        for group in groups:
            group = list(group) #List of tuples
            #Check if group is a cutset
            aux_graph = g.copy()
            aux_graph.delete_edges(group)
            cut_set = True
            for node in nodes:
                cut_set = cut_set and aux_graph.vertex_disjoint_paths(s, node, neighbors = "ignore") == 0 #There is no path between s and node               
            if cut_set:
                return i
    return e


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

def removeLinks(graph, vertex):
    """
    return: A new graph where vertex is disconnected
    """
    new_graph = graph.copy()
    neighbors = new_graph.neighbors(vertex)
    for target in neighbors:
        new_graph.delete_edges([(vertex, target)])
    return new_graph

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

def sortEigenValuesVectors(eigenvalues, eigenvectors, desc=False,absValue= True):
    """
    Auxiliary function
    eigenvalues: Array of eigenvalues
    eigenvectors: Array of eigenvector
    asc: Sort reverse
    abs: Take tha absolute value to compare
    return: Sorted eigen values and eigen vectors, eigenvectors[:,i] is the eigenvector corresponding to the eigenvalues[i]
    """

    map_values = {}
    for value in eigenvalues:       
        map_values[abs(value)] = value
    if absValue:
        eigenvalues = np.abs(eigenvalues)  
    eigenvalues = list(eigenvalues)
    eigenvectors = list(eigenvectors)

    pairs = list(zip(eigenvalues, eigenvectors))
    pairs = sorted(pairs, reverse=desc, key=lambda x: x[0])
    unzip = list(zip(*pairs))
    eigenvalues = np.array(unzip[0])
    eigenvectors = np.array(unzip[1])
    real_eigenvalues = []
    for value in eigenvalues:
        real_eigenvalues.append(map_values[value])
    return real_eigenvalues, eigenvectors

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

def makeRandomServices(graph):
    """
    retrun: random set of services available locally and the set of nonlocal services for each vertex. List = [[[A(v_0)], [N(v_0)]], ... , [[A(v_n-1)], [N(v_n-1)]]]
    """
    services = []
    v = graph.vcount()
    vertices = list(range(v))
    for i in range(v):
        non_local = []
        amount_local = random.randint(0, v-1)
        local = random.sample(vertices, amount_local)
        available_services = [service for service in vertices if service not in local]
        amount_nonlocal = random.randint(1, v - amount_local)
        non_local = random.sample(available_services, amount_nonlocal)
        services.append([local, non_local])
    return services

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

def generateWeight(g, edge =True, vertex=False, name= 'weight', negative = False):
    """
    g: Graph
    negative: Set negative values
    returns: Graph with random weights on its edges, vertex or both
    """
    if edge:
        w = np.random.rand(g.ecount())
        if negative:
            w = -w
        g.es[name] = w
    if vertex: 
        w = np.random.rand(g.vcount())
        if negative:
            w = -w
        g.vs[name] = w
    return g

def generatePositions(g, x_position='x_position', y_position='y_position'):
    """
    returns graph with vertex attribute x_position and y_position
    """
    n = g.vcount()
    x = np.random.rand(n)
    y = np.random.rand(n)
    g.vs[x_position] = x
    g.vs[y_position] = y
    return g

def generateDistance(g, x_position, y_position, distance):
    """
    distance: Name of distance attribute
    returns graph with edge attribute that indicates the distance of the nodes
    """
    distances = []
    for edge in g.es:
        x1 = g.vs[x_position][edge.source]
        y1 = g.vs[y_position][edge.source]
        x2 = g.vs[x_position][edge.target]
        y2 = g.vs[y_position][edge.target]
        d = math.sqrt(((x1 - x2)**2 + (y1 - y2)**2))
        distances.append(d)
    g.es[distance] = distances
    return g

def setStates(g):
    """
    g: Graph
    returns: Graph with random states between [0, 1] for each vertex in attribute 'state'
    """
    s = np.random.uniform(0,1, g.vcount())
    s[s > 0.999] = 1
    g.vs['state'] = s
    return g

def normalize(v):
    """
    v: List of eigenvectors
    return: Normalized eigenvectors
    """
    n = np.linalg.norm(v, axis=1)
    return v/n

def getDAG(g, reverse = False):
    """
    g: Graph
    Return directed acyclic graph
    """
    n_graph = nx.DiGraph()
    n_graph.add_nodes_from(list(range(g.vcount())))
    for edge in g.es:
        if reverse:
            target = edge.source
            source = edge.target
        else:
            source = edge.source
            target = edge.target
        n_graph.add_edge(source, target)
        try:
            nx.find_cycle(n_graph)
            #There is a cycle
            n_graph.remove_edge(source, target)
        except:
            pass
    graph = Graph(directed=True)
    graph.add_vertices(g.vcount())
    edges = list(n_graph.edges)
    graph.add_edges(edges)

    return graph

def getPaths(g, M, u):
    """
    g: Graph
    M: List of vertex
    u: List of vertex
    return: List of all paths where source vertex belongs to M and target vertex belongs to u
    """
    paths = []
    n_graph = nx.DiGraph()
    n_graph.add_nodes_from(list(range(g.vcount())))
    n_graph.add_edges_from(g.get_edgelist())

    for i in M:
        for j in u:
            paths_ij = list(nx.all_simple_paths(n_graph,source=i,target=j))
            for path in paths_ij:
                paths.append(path)

    return paths

def h_vi(v, paths):
    """
    v: Vertex
    paths: List of paths
    """
    k= 0
    M = len(paths)
    for p in paths:
        if p[0] == v:
            k+= 1
    return k * math.log(1/M) / M

def H_f(g):
    """
    g: Graph
    return: Measure for treeness
    """
    k_i = np.array(g.indegree())
    k_out = np.array(g.outdegree())
    M = np.where(k_i == 0)[0]
    u = np.where(k_out == 0)[0]
    acc = 0
    paths = getPaths(g, M, u)

    for vi in M:
        acc += h_vi(vi, paths)

    return acc / len(M)

def criticality(graph, n):
    """
    n: Iteration number
    return: Criticality of all vertices at iteration n
    """
    h_n_1 = np.ones(graph.vcount())
    h_n = np.zeros(graph.vcount())
    for it in range(n-1):
        for v in range(graph.vcount()):
            acc = 0
            v_neighbors = graph.neighbors(v)
            k_j = len(v_neighbors)
            for neighbor in v_neighbors:
                acc += h_n_1[neighbor] / k_j
            h_n[v] = acc
        h_n_1 = h_n
        h_n = np.zeros(graph.vcount())
    return h_n_1

def rth_nearest_neighbors(graph, r):
    """
    return: the number of rth nearest neighbors of all nodes in graph
    """
    n = graph.vcount()
    reachable = graph.neighborhood(order=r) #array
    result = np.zeros(n)
    for vertex in range(n):
        for neighbor in reachable[vertex]:
            if graph.shortest_paths_dijkstra(vertex, neighbor)[0][0] == r:
                result[vertex] += 1
    return result

def average_rth_nearest_neighbors(graph, r):
    """
    return: average number of rth-nearest neighbors of all nodes in the graph
    """
    return np.sum(rth_nearest_neighbors(graph, r)) / graph.vcount()
            
def getGiantComponent(graph):
    """
    return: igraph object of the giant component of graph
    """
    components = graph.components()
    giant_component_list = []
    size = 0
    for component in components:
        if len(component) > size:
            size = len(component)
            giant_component_list = component
    
    #Add edges to giant component
    giant_component = Graph()
    giant_component.add_vertices(len(giant_component_list))
    for i in range(len(giant_component_list)): #i new vertex id
        neighbors = graph.neighbors(giant_component_list[i])
        for neighbor in neighbors:
            if neighbor in giant_component_list:
                new_index = giant_component_list.index(neighbor) #new index of target
                if giant_component.get_eid(i, new_index, directed=False, error=False) == -1:
                    giant_component.add_edge(i, new_index)
    return giant_component

def localConnectivityAux(graph):
    rm = 0
    n = graph.vcount()
    for i in range(n):
        for j in range(n):
            if i != j:
                shortest_path_length = graph.shortest_paths_dijkstra(i, j)[0][0]
                if shortest_path_length > rm and not np.isinf(shortest_path_length):
                    rm = shortest_path_length
    
    acc = 0
    for i in range(1, rm + 1):
        acc += average_rth_nearest_neighbors(graph, i) / i
    return acc

def performance(graph):
    """
    Auxiliary function for vulnerability
    returns the performance of a graph
    """
    acc = 0
    n = graph.vcount()
    for i in range(n):
        for j in range(n):
            if i != j:
                shortest_path_length = graph.shortest_paths_dijkstra(i, j)[0][0]
                if not np.isinf(shortest_path_length):
                    acc += shortest_path_length
    try:
        return acc / (n * (n-1))
    except ZeroDivisionError:
        return None

def attack_edges(graph, n= None):
    """
    returns a new graph with n random edges deleted, if n not specified, then n is chosen randomly
    """
    if n is None:
        n_edges = graph.ecount()
        n = random.randint(1, n_edges)
    attacked_graph = graph.copy()
    for i in range(n):
        edges = list(range(attacked_graph.ecount()))
        edge_to_delete = random.choice(edges)
        attacked_graph.delete_edges(edge_to_delete)   
    return attacked_graph

def setRoadState(graph, p):
    """
    Set edge attribute edges with attribute 'state', state can be open or closed. States are chosen randomly
    p: Probability of closed road
    """
    states = []
    for i in range(graph.ecount()):
        n = random.random()
        if n <= p:
            states.append('closed')
        else:
            states.append('open')
    graph.es['state'] = states
    return graph

def generateVoltages(g, name, generation_nodes):
    """
    generation_nodes: List of generation nodes, the voltage of this nodes is positive
    returns a graph with attribute name on nodes, where each value is the voltage of the node
    """
    n = g.vcount()
    w = -np.random.rand(n)
    for gen_node in generation_nodes:
        w[gen_node]  = -w[gen_node]
    g.vs[name] = w
    return g

def setId(g, name='id', vertex = True, edge=False):
    if vertex:
        node_id = list(range(g.vcount()))
        g.vs[name] = node_id
    if edge:
        edge_id = list(range(g.ecount()))
        g.es[name] = edge_id
    return g


def communicationNetworkLoad(g):
    """
    Returns an array with the load of each node in the communication network g
    """
    n = g.vcount()
    load = np.zeros(n)    
    for i in range(n):
        sum_i = 0
        for j in range(n):
            for k in range(n):
                if i != j and j != k and i != k:
                    n_i_jk = 0
                    all_shortest_paths = g.get_all_shortest_paths(j, k)
                    n_jk = len(all_shortest_paths)
                    for path in all_shortest_paths:
                        if i in path:
                            n_i_jk += 1
                    sum_i += n_i_jk / n_jk
        load[i] = sum_i
    return load

def powerGridNetworkLoad(g, voltage, admittance):
    """
    voltage: Node attribute for voltage
    admittance: Edge attribute for admittance
    Returns an array with the load of each node
    """
    n = g.vcount()
    adjacency_matrix = -np.array(g.get_adjacency(attribute=admittance).data)
    #Get voltages
    voltages = np.array(g.vs[voltage])
    #Change adjacency
    for i in range(n):
        if voltages[i] > 0:
            row = np.zeros(n)
            row[i] = 1
            adjacency_matrix[i] = row
        else:
            adjacency_matrix[i][i] = -np.sum(adjacency_matrix[i])
    current_in_nodes = np.matmul(adjacency_matrix, voltages)
    nodes_power_load = voltages * current_in_nodes
    return nodes_power_load

def getGiantComponentList(graph):
    """
    return: List of nodes in the giant connected component
    """
    components = graph.components()
    giant_component_list = []
    size = 0
    for component in components:
        if len(component) > size:
            size = len(component)
            giant_component_list = component
    return giant_component_list

def cascadingFailures(CN, PG, max_loads_CN, max_loads_PG, voltage, admittance, CN_dict, PG_dict):
    """
    returns: The number of failed nodes
    """
    n = 0
    big_loop = True
    while(big_loop):
        big_loop = False
        sub_graph_detection = True
        while(sub_graph_detection):
            sub_graph_detection = False
            #Check CN
            CN_components = CN.components()
            if len(CN_components) > 1:
                sub_graph_detection = True
                giant_component = getGiantComponentList(CN)
                n_CN = CN.vcount()
                #Delete nodes outside the giant component
                for v in reversed(range(n_CN)):
                    if not v in giant_component:
                        n += 1
                        original_id = CN.vs['id'][v]
                        CN.delete_vertices(v)
                        #Delete their counterpart nodes in the PG
                        if original_id in CN_dict:
                            n += 1
                            PG_neighbor_id = CN_dict[original_id]
                            PG_nodes_id = PG.vs['id']
                            PG_neighbor = PG_nodes_id.index(PG_neighbor_id)
                            #Delete from graph
                            PG.delete_vertices(PG_neighbor)
                            #Delete nodes from dictionaries
                            del CN_dict[original_id]
                            del PG_dict[PG_neighbor_id]

            #Check PG
            PG_components = PG.components()
            components_to_delete = []
            for component in PG_components:
                delete_component = True
                for v in component:
                    delete_component = delete_component and PG.vs[v][voltage] < 0
                if delete_component:
                    components_to_delete += component
            if len(components_to_delete) > 0:
                sub_graph_detection = True
                n += len(components_to_delete)
                #Delete their counterpart in CN
                for v in components_to_delete:
                    original_id = PG.vs['id'][v]
                    if original_id in PG_dict:
                        n += 1
                        CN_neighbor_id = PG_dict[original_id]
                        CN_nodes_id = CN.vs['id']
                        CN_neighbor = CN_nodes_id.index(CN_neighbor_id)
                        CN.delete_vertices(CN_neighbor)
                        del CN_dict[CN_neighbor_id]
                        del PG_dict[original_id]
                PG.delete_vertices(components_to_delete)               

        #Load redistribution
        #CN load
        CN_load = communicationNetworkLoad(CN)
        CN_nodes_to_delete = []
        for i in range(len(CN_load)):
            actual_load = CN_load[i]
            original_id = CN.vs[i]['id']
            if actual_load > max_loads_CN[original_id]:
                CN_nodes_to_delete.append(i)
        if len(CN_nodes_to_delete) > 0:
            big_loop = True
            n += len(CN_nodes_to_delete)
            #Delete their counterpart in PG
            for v in CN_nodes_to_delete:
                original_id = CN.vs['id'][v]
                if original_id in CN_dict:
                    n += 1
                    PG_neighbor_id = CN_dict[original_id]
                    PG_nodes_id = PG.vs['id']
                    PG_neighbor = PG_nodes_id.index(PG_neighbor_id)
                    PG.delete_vertices(PG_neighbor)
                    del PG_dict[PG_neighbor_id]
                    del CN_dict[original_id]
            CN.delete_vertices(CN_nodes_to_delete)           

        #PG load
        PG_load = powerGridNetworkLoad(PG, voltage, admittance)
        PG_nodes_to_delete = []
        for i in range(len(PG_load)):
            actual_load = PG_load[i]
            original_id = PG.vs[i]['id']
            if actual_load > max_loads_PG[original_id]:
                PG_nodes_to_delete.append(i)
        if len(PG_nodes_to_delete) > 0:
            big_loop = True
            n += len(PG_nodes_to_delete)
            #Delete their counterpart in CN
            for v in PG_nodes_to_delete:
                original_id = PG.vs['id'][v]
                if original_id in PG_dict:
                    n += 1
                    CN_neighbor_id = PG_dict[original_id]
                    CN_nodes_id = CN.vs['id']
                    CN_neighbor = CN_nodes_id.index(CN_neighbor_id)
                    CN.delete_vertices(CN_neighbor)
                    del CN_dict[CN_neighbor_id]
                    del PG_dict[original_id]
            PG.delete_vertices(PG_nodes_to_delete)
                
    return n

def toDict(interdependencies):
    """
    interdependencies: List of interdependencies, format: [(CN, PG)]
    Important: Each node can have at most one interdependent link
    return: two dictionaaries, one for each network, that indicates the neighbors in the other network
    """
    CN_dict = {}
    PG_dict = {}
    for edge in interdependencies:
        CN_v = edge[0]
        PG_v = edge[1]
        #Add to CN_dict
        CN_dict[CN_v] = PG_v
        #Add to PG_dict
        PG_dict[PG_v] = CN_v
    return CN_dict, PG_dict

def spanningTrees(g):
    """
    return: List of all spanning trees of graph g, each list represent edges id
    """
    m = g.ecount()
    n = g.vcount()
    edges = list(range(m))
    comb = combinations(edges, n -1)
    spanning_trees = []
    for group in comb:
        group = list(group)
        # Make graph
        edge_list = []
        for edge in group:
            source = g.es[edge].source
            target = g.es[edge].target
            edge_list.append((source, target))
        aux_g = nx.from_edgelist(edge_list)
        # Check cycle
        try:
            nx.find_cycle(aux_g)
        except:
            spanning_trees.append(group)
    return spanning_trees

def cospanningTrees(g):
    """
    return: List of all cospanning trees of graph g, each list represent edges id
    """
    spanning_trees = spanningTrees(g)
    cospanning_trees = []
    m = g.ecount()
    for st in spanning_trees:
        ct = []
        for e in range(m):
            if e not in st:
                ct.append(e)
        cospanning_trees.append(ct)
    return cospanning_trees

