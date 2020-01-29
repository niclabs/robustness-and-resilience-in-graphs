from distanceBasedCharacteristics import deltaEfficiency
from igraph import *

g = Graph.Erdos_Renyi(n=10, m=20)
print(deltaEfficiency(g))