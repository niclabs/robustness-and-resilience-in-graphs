import numpy as np
from igraph import *
import random
import itertools
import math
from auxiliaryFunctions import *
import warnings

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

def preferentialPerturbation(g1, g2):
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