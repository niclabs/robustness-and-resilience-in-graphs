import numpy as np
from igraph import *
import random
import itertools
import math
from auxiliaryFunctions import *
import warnings

def normalizedGCC(graph, t=0, state=None, distance=None, resourceCenter=0, p=0.3):
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

def splittingNumber(g, k=2, numOfLoops=1, seed=1, dev= 1):
    """
    g: Graph
    k: Number of components   
    numOfLoops: Parameter for the algorithm, default = 1
    seed: Seed for randomize, if seed = 0, this parameter is not used, default=0
    dev: Parameter for the algorithm, standard deviation
    return: The average number of edges that need to be removed to break  g into k connected components, -1 if it's not possible
    """
    givenComponents = k
    i = 0
    cond = givenComponents >= len(g.components())

    l = g.ecount()
    s = 0
    for j in range(1, l + 1):
        s += math.factorial(l) / (math.factorial(i) * math.factorial(l - i))
    numOfLinks = np.zeros(int(s) + 1)

    if(seed):
        random.seed(seed)

    while (cond):
        auxGraph = g.copy() #make all links functioning
        numOfComponents = len(auxGraph.components()) #Calculate numOfcomponents

        while (numOfComponents < givenComponents):
            #Choose randomly an uninterrupted link
            numberTotalLinks = auxGraph.ecount()
            linkId = random.randint(0, numberTotalLinks - 1)            

            numOfLinks[i] += 1
            auxGraph.delete_edges([linkId]) #interrup the link
            numOfComponents = len(auxGraph.components())#compute numOfComponents

            if(i >= 1 and i >= numOfLoops):
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    mean1 = np.mean(numOfLinks[1:i+1]) #mean value of numOfLinnks[1:i]
                    mean2 = np.mean(numOfLinks[1:i]) #mean value of numOfLinnks[1:i - 1]
                    if(abs(mean1 - mean2) < dev):
                        return (mean1 + mean2)/2
        i += 1
    return -1

def randomRobustnessIndex(g, m=1):
    """
    g: Graph
    m: number of rounds of attacks, default=1
    """
    v = g.vcount()
    acc = 0
    for i in range(m): #Each random attacks
        for j in range(1, v+1): #Remove j random vertices
            aux = g.copy()
            vertices = list(range(aux.vcount()))
            vertices_to_delete = random.sample(vertices, j)
            aux.delete_vertices(vertices_to_delete)
            clusters = aux.components()
            maxCluster = 0
            for c in clusters:
                l = len(c)
                if l > maxCluster:
                    maxCluster = l
            acc += maxCluster / v
    try:
        result =  (1/m) * (1/v) * acc
    except ZeroDivisionError:
        result = None
    return result

def robustnessMeasure53(g):
    """
    g: Graph
    return: Messure 5.3
    """
    aux = g.copy() #The function doesn't modify the graph
    n = aux.vcount()
    if (n < 1):
        return -1 #Error
    Ct = 0
    for i in range(n):
        degrees = aux.degree()
        v = np.argmax(degrees) #Vertex of maximum degree
        aux.delete_vertices([v])
        C = aux.components() #Components of the graph
        Ci = 0 #Order of the biggest component
        for comp in C:
            compOrder = len(comp)
            if(compOrder > Ci):
                Ci = compOrder 
        Ct += Ci
    try: 
        result = Ct / n
    except ZeroDivisionError:
        result = None
    return result

def connectivityRobustnessFunction(g, k=1):
    """
    g: Graph
    k: Number of vertices removed, default = 1
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

def kResilienceFactor(g, k=2):
    """
    g: Graph
    k: k - 1 vertices will be removed, default k=2
    return: The percentage of connected components of g that remain connected after the removal k - 1 vertices 
    """
    n = g.vcount()
    if(k - 1 > n):
        return None #Error
       
    vertices_list = np.arange(n) #List of vertices
    combinations = list(itertools.combinations(vertices_list, k -1))
    connected = float(0)

    for vertices in combinations: #Remove vertices
        aux_graph = g.copy()
        aux_graph.delete_vertices(vertices)
        if len(aux_graph.components()) == 1:
            connected += 1
    return connected / len(combinations)

def resilienceFactor(g):
    n = g.vcount()
    if( n < 3):
        return None #Error
    result = np.zeros(n-2)
    for i in range(2, n):
        auxGraph = g.copy()
        result[i-2] = kResilienceFactor(auxGraph, k=i)
    return np.mean(result)

def perturbationScore(g, p=perturbationFunction):
    """
    g: Graph
    p: Perturbation function, returns a graph
    """
    aux = p(g)
    return perturbationScoreTwo(g, aux)

def preferentialPerturbation(g1, g2): # auxiliary
    """
    g1: Graph
    g2: Graph
    return: A list of vertices that has to be removed in order to maximize the pertrbation in g1 and minimize the perturbation in g2
    """
    max = min(g1.vcount(), g2.vcount()) #Maximum amount of vertices that  can be removed
    nodes = [] #Removed vertices in the perturbation
    perturbation1 = float('-inf')
    perturbation2 = float('inf')
    nodes = list(range(max))
    for i in range(1, max + 1):
        comb = list(itertools.permutations(nodes, i))
        for l_vertex in comb: #l_vertex: tuples
            auxg1 = g1.copy()
            auxg2 = g2.copy()
            auxg1.delete_vertices(l_vertex)
            auxg2.delete_vertices(l_vertex)
            partial_p1 = perturbationScoreTwo(g1, auxg1)
            partial_p2 = perturbationScoreTwo(g2, auxg2)
            if(partial_p1 > perturbation1 and  partial_p2 < perturbation2):
                nodes = np.array(l_vertex)
                perturbation1 = partial_p1
                perturbation2 = partial_p2
    return nodes


def maximumPerturbationScore(g1, g2):
    """
    g1: Graph
    g2: Graph
    return: The ratio between the preferential perturbation score and the proportion of vertices that are removed,
            it returns -1 when vertices are not removed
    """
    p = preferentialPerturbation(g1, g2)
    n = len(p)
    if n == 0:
        return -1
    auxg1 = g1.copy()
    auxg2= g2.copy()
    auxg1.delete_vertices(p)
    auxg2.delete_vertices(p)
    p1 = perturbationScoreTwo(g1, auxg1)
    p2 = perturbationScoreTwo(g2, auxg2)
    p12 = p1 - p2
    n1 = g1.vcount()
    return p12 / (n / n1)

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
    vertex_id = list(range(n))
    for vertex in ranking:
        # Delete vertex
        id = vertex_id.index(vertex)
        graph.delete_vertices(id)
        vertex_id.remove(vertex)
        # Calculate the size of giant connected component
        comp = graph.components()
        sizes = comp.sizes()
        if len(sizes) != 0:
            acc += max(sizes)
    try:
        result = acc / n
    except ZeroDivisionError:
        result = None
    return result
