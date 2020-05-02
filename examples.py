from igraph import *
from otherCharacteristics import entropy

g = Graph.Erdos_Renyi(n=20, m=30)
print(len(g.components()))
print(entropy(g))
