import numpy as np
from igraph import *
import random

def splittingNumber(g, k, dev= 1):
    """
    g: Graph
    k:
    dev:
    return: The average number of edges that need to be removed to break  g into k connected components
    """

    numOfLoops = 0 #TODO: definir
    givenComponents = k

    i = 0
    numOfComponents = 0
    cond = True
    l = g.ecount()
    k = 3 #TODO: Calculate
    numOfLinks = np.zeros(k)
    while (cond):
        auxGraph = g.copy() #make all links functioning
        numOfComponents = len(auxGraph.components()) #Calculate numOfcomponents
        numOfLinks[i] = 0

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

def randomRobustnessIndex(g, m):
    """
    g: Graph
    m:
    return
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
    return
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

def connectivityRobustnessFunction(g, k):
    """
    g: Graph
    k: Number of vertices removed
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

def kResilienceFactor(g, k):
    """
    g: Graph
    k:
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
        result[i-2] = kResilienceFactor(auxGraph, i)
    return np.mean(result)

def pertubationScore(g, p):
    """
    g: Graph
    p: Perturbation function
    return:
    """
    BComp = g.components()
    before = 0
    for c in BComp:
        size = len(c)
        if(size > before):
            before = size
    aux = p(g)
    AComp = aux.components()
    after = 0
    for c in AComp:
        size = len(c)
        if(size > after):
            after = size
    return (before - after)/ before

def preferentialPerturbation(g1, g2):
    """
    g1: Graph
    g2: Graph
    return:
    """
    return

def maximunPerturbationScore(g1,g2):
    """
    g1: Graph
    g2: Graph
    return:
    """
    return