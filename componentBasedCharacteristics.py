import numpy as np
from igraph import *
import random
import itertools
import math

def splittingNumber(g, k, numOfLoops=1, seed=0, dev= 1):
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
    cond = True

    l = g.ecount()
    s = 0
    for j in range(1, l + 1):
        s += math.factorial(l) / (math.factorial(i) * math.factorial(l - i))
    numOfLinks = np.zeros(s)

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

            if(i > 1 and i > numOfLoops):
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
    return: Messure 5.3
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
        result[i-2] = kResilienceFactor(auxGraph, k=i)
    return np.mean(result)

def perturbationScore(g, p):
    """
    g: Graph
    p: Perturbation function, returns a graph
    """
    aux = p(g)
    return perturbationScoreTwo(g, aux)

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


def maximunPerturbationScore(g1, g2):
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