from igraph import *
from characteristicsForRankingElements import *
g = Graph.Erdos_Renyi(n=20, m=30)

#g = Graph([(0,1), (1,2), (2,3), (3,0)])
print(mcv(g))


