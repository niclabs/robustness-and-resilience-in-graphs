from igraph import *
from characteristicsForRankingElements import *
g = Graph([(0,2), (2,4), (1,4), (3,4)])
#g = Graph([(0,1)])
print(mcv2(g))