from igraph import *
import numpy as np
from spectralCharacteristics import *
g = Graph.Erdos_Renyi(n=200, m=210)
#g = Graph([(0,3), (1,2), (1,3), (2,4), (3,6), (6,5), (3,5)])
print(generalizedRobustnessIndex2(g, 30))
