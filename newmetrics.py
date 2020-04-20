from igraph import *
import numpy as np
from auxiliaryFunctions import *

def percentageOfNoncriticalNodes(CN, PG, interdep= None, p=0.85, admittance = None, voltage = None, tolerance = 0.5, a_t = 0.2, a_p = 0.2):
    """
    Note: CN and PG must be the same number of nodes
    CN: Communication network
    PG: Power grid network
    interdep: List of interdependencies, format [(CN, PG)], each node can have at most one interdependent link
    p: Percentage of interdependent nodes, only used when interdep is None
    admittance: Edge attribute for admittance
    voltage: Vertex attribute for voltage
    generation_nodes: List of generation nodes
    tolerance: Tolerance threshold
    a_t: Tolerance parameter of communication network 
    a_p: Tolerance parameter of nodes in power grid network

    """
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

def robustnessIndex(g, rank=None):
    """
    g: Graph
    rank: Rank function, it gives a list of nodes in order of importance, default: Pagerank
    """
    rank_list = []
    n = g.vcount()
    if rank is None:
        rank_list = np.array(g.pagerank())
        rank_list = rank_list.argsort()[::-1][:n]
    else:
        rank_list = rank(g)
    
    g_aux = g.copy()
    # Remove nodes and calculate the size of giant connected component
    for i in range(n - 1):
        vertex = rank_list[i]
        neighbors = g_aux.neighbors(vertex, mode='OUT')
        # Delete edges from vertex
        for target in neighbors:
            g_aux.delete_edges([(vertex, target)])
        # Get size of giant component
        comp = g_aux.components()
        rank_list[i] = max(comp.sizes())
    rank_list[n-1] = 0

    return np.sum(rank_list) / n

def treeness(g):
    """
    g: Graph
    """
    Gc = getDAG(g)
    Gc_r = getDAG(g, reverse=True)
    hf = H_f(Gc)
    hb = H_f(Gc_r)
    return hf- hb / max(hf, hb)

def entropy(graph, gamma= 0.5, x_pos=None, y_pos=None):
    """
    gamma: Weight of the node criticality, weight of edge criticality 1 - gamma
    x_pos: Parameter that indicates the name attribute for x position of vertices
    y_pos: Parameter that indicates the name attribute for y position of vertices
    If x_pos and y_pos are not indicated will be set randomly and returned with the final result
    """
    if x_pos is None or y_pos is None:
        x_pos = np.random.rand(graph.vcount())
        y_pos = np.random.rand(graph.vcount())    
    else:
        x_pos = graph.vs[x_pos]
        y_pos = graph.vs[y_pos]
    # Calculate node criticality
    CV = criticality(graph, graph.vcount())
   
    edge_correlation_factor = np.zeros((graph.vcount(), graph.vcount()))
    # Calculate edge correlation factor matrix
    for i in range(graph.vcount()):
        for j in range(graph.vcount()):
            if i == j:
                edge_correlation_factor[i][j] = 1
            else:
                edge_correlation_factor[i][j] = max(CV[i], CV[j]) / (CV[i] + CV[j])
    
    adjacency = graph.get_adjacency().data
    # Change adjacency matrix
    for i in range(graph.vcount()):
        for j in range(graph.vcount()):
            if i == j:
                adjacency[i][j] = 1
    
    average_transmission_efficiency = np.zeros(graph.vcount())
    # Calculate average transmission efficiency vector
    for i in range(graph.vcount()):
        acc = 0
        for j in range(graph.vcount()):
            if i!= j:
                distance = 0 #TODO: Calculate
                acc += 1 / distance
        average_transmission_efficiency[i] = (2 * acc) / (graph.vcount() * (graph.vcount() - 1))
    
    # Average transmission efficiency matrix
    I = np.zeros((graph.vcount(), graph.vcount()))
    for i in range(graph.vcount()):
        for j in range(graph.vcount()):
            I[i][j] = average_transmission_efficiency[i]
    
    # Calculate edge criticality matrix
    CE = I * adjacency * edge_correlation_factor

    # Calculate comprehensive critical degree of nodes
    CS = np.zeros(graph.vcount())
    for i in range(graph.vcount()):
        acc = 0
        i_neighbors = graph.neighbors(i)
        for j in i_neighbors:
            acc += CE[i][j]
        CS[i] = (gamma * CV[i]) + (1 - gamma) * acc / len(i_neighbors)

    # Calculate critical coefficient of nodes
    S = CS / np.sum(CS)
    result = - np.sum( S * np.log(S))
    return result, x_pos, y_pos

def globalResilienceAnalysis():
    pass

def robustnessMeasureR(graph, ranking_function=None):
    """
    ranking_function: Function that returns a list of ranked nodes in graph, use: ranking_function(graph), if ranking_function is None, then the nodes are ranked by degree
    """
    if ranking_function is None:
        degrees = graph.degree()
        ranking = np.flip(np.argsort(degrees))
    else:
        ranking = ranking_function(graph)
    acc = 0
    n = graph.vcount()
    for vertex in ranking:
        # Delete vertex
        graph.delete_vertices(vertex)
        # Calculate the size of giant connected component
        comp = graph.components()
        acc += max(comp.sizes())
    try:
        result = acc / n
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

def localConnectivity(graph):
    giant_component = getGiantComponent(graph)
    try:
        result = localConnectivityAux(graph) / localConnectivityAux(giant_component)
    except ZeroDivisionError:
        result = None
    return result

def vulnerability(graph, attack_function=attack_edges):
    """
    attack_function: Function that damages graph by deleting vertices, edges or both, returns a graph, use: attack_function(graph)
    """
    performance_normal = performance(graph)
    attacked_graph = attack_function(graph)
    performance_damage = performance(attacked_graph)
    try:
        result = (performance_normal - performance_damage) / performance_damage
    except ZeroDivisionError:
        result = None
    return result

def normalizedGiantConnectedComponent(graph, t=0, state=None, distance=None, resourceCenter=0, p=0.3):
    """
    t: Time, starts at 0 where all closed/disrupted edges are removed
    state: Edge attribute that indicates if an edge is open or closed   
    distance: Edge attribute that indicates the distance of each road
    resourceCenter: Node id that indicates the resource center node
    p: Probability that a road is closed, used when states are not asigned
    """
    # Check n >0 and m >0
    if graph.vcount() == 0 or graph.ecount() == 0:
        return None
    comp = graph.components()
    gcc_baseline = max(comp.sizes())

    # Set state
    if state is None:
        state = 'state'
        graph = setRoadState(graph, p)
    
    # Set distance
    if distance is None:
        distance = 'distance'
        graph = generateWeight(graph, edge=True, vertex=False, name=distance)
    edges_to_add = []
    # Identify closed/disrupted edges
    for e in graph.es:
        if e[state] != 'open':
            edge = (e.source, e.target)
            edges_to_add.append(edge)
    # There is no edge to delete
    if t >= len(edges_to_add):
        return 1.0
    
    else:
        # Rank closed/disrupted edges
        rank_graph = graph.copy()
        rank_graph.delete_edges(edges_to_add)
        shortestDistanceToEdge = []
        for pair in edges_to_add:
            v1 = pair[0]
            v2 = pair[1]
            d1 = rank_graph.shortest_paths_dijkstra(resourceCenter, v1, weights=distance)[0][0]
            d2 = rank_graph.shortest_paths_dijkstra(resourceCenter, v2, weights=distance)[0][0]
            shortestDistanceToEdge.append(min(d1, d2))
        argsorted_distances = np.argsort(shortestDistanceToEdge)
        sorted_edges_to_add = list()
        for i in argsorted_distances:
            sorted_edges_to_add.append(edges_to_add[i])

        # delete closed/disrupted edges
        edges_to_delete = sorted_edges_to_add[t:]
        aux_graph = graph.copy()
        aux_graph.delete_edges(edges_to_delete)
        comp_aux = aux_graph.components()
        gcc_t = max(comp_aux.sizes())
        return gcc_t / gcc_baseline