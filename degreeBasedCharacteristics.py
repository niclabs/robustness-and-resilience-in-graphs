import math
import numpy as np
from igraph import *
import sys
from auxiliaryFunctions import *

def degreeEntropy(g):
    """
    g: Graph
    """
    sum = 0
    for i in range(1, g.vcount()):
        p_i = getProbabilityDegree(g, i)
        if(p_i != 0):
            sum += p_i * math.log(p_i)
    return -sum

def relativeEntropy(g):
    """
    g: Graph
    """
    n = g.vcount()
    pk = getDegreeDistribution(g)
    sum = 0
    for i in pk:
        if (n * i == 0):
            sum += i * math.log(sys.float_info.min)
        else:
            sum += i * math.log(n * i)
    return sum