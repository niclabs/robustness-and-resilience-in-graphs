from igraph import *
import numpy as np
from pathBasedCharacteristics import *
#g = Graph.Erdos_Renyi(n=200, m=210)

g = Graph.Barabasi(15,3)
print(g.vcount())
l = makeEmptyServices(g)
print(len(l))
print(l)
