from igraph import *
import numpy as np
from pathBasedCharacteristics import *
#g = Graph.Erdos_Renyi(n=200, m=210)

g = Graph.Barabasi(15,3)


#Ejemplo de paper
g = Graph([(0,2), (0,1), (0,3), (1,3), (1,2), (2,3), (2,4), (4,5), (5,6), (5,7), (4,6), (4,7), (6,7)])
A0  = [[0,1,2], [3]]
A1 = [[0,1,3], [2]]
A2 = [[0,2,3], [1]]
A3 = [[1,2,3], [0]]
A4 = [[4,5,6], [7]]
A5 = [[4,5,7], [6]]
A6 = [[4,6,7], [5]]
A7 = [[5,6,7], [4]]

L = [A0, A1, A2, A3, A4, A5, A6, A7]

print(vertexResilience(g,L))
