from igraph import *
import unittest
import numpy as np
from characteristicsForRankingElements import *

class VertexLoadTest(unittest.TestCase):
    def setUp(self):
        self.directed = Graph([(1,2), (1,3), (3,4)], directed =True)
        self.g =  Graph([(1,2), (1,3), (3,4)])
    
    def testDirectedGraph(self):
        self.assertEqual(vertexLoad(self.directed, v=0, n=1), 0, "Error in diredted graph, case v has no neighbors")
        self.assertEqual(vertexLoad(self.directed, v=0, n=0), 1, "Error in diredted graph, case v has no neighbors and n= 0")
        self.assertEqual(vertexLoad(self.directed, v=0, n=5), 0, "Error in diredted graph, case v has no neighbors and n= 5")
        self.assertEqual(vertexLoad(self.directed, v=1, n=1), 6, "Error in directed graph, case v has neighbors anmd n= 1")
        self.assertEqual(vertexLoad(self.directed, v=1, n=0), 1, "Error in directed graph, case v has neighbors and n= 0")
        self.assertEqual(vertexLoad(self.directed, v=1, n=5), 7776, "Error in directed graph, case v has neighbors and n= 6")

    def testGraph(self):
        self.assertEqual(vertexLoad(self.g, v=0, n=1), 0, "Error in graph, case v has no neighbors")
        self.assertEqual(vertexLoad(self.g, v=0, n=0), 1, "Error in graph, case v has no neighbors and n= 0")
        self.assertEqual(vertexLoad(self.g, v=0, n=5), 0, "Error in graph, case v has no neighbors and n= 5")
        self.assertEqual(vertexLoad(self.g, v=1, n=1), 6, "Error in graph, case v has neighbors anmd n= 1")
        self.assertEqual(vertexLoad(self.g, v=1, n=0), 1, "Error in graph, case v has neighbors and n= 0")
        self.assertEqual(vertexLoad(self.g, v=1, n=5), 7776, "Error in graph, case v has neighbors and n= 6")

class RandomWalkTest(unittest.TestCase):
    def setUp(self):
        self.directed = Graph([(0,1), (0,6), (1,4), (1,3), (2,3), (3,1), (1,2), (5,1)], directed=True)
        self.directedOneWay = [0, 6]
        self.directedTwoWays = [0, 1, 3]
        self.seed = 1

        self.g = Graph([(1,2), (5,4), (4,7), (4,6), (5,6), (6,7)])
        self.oneWay = [1,2]
        self.twoWays = [5, 4, 7]
        

    def testDirectedGraph(self):
        oneWay = randomWalk(self.directed, 0, 6)
        self.assertEqual(oneWay, self.directedOneWay)
        empty = randomWalk(self.directed, 0, 5)
        self.assertEqual(empty, [])
        twoPossibleWays = randomWalk(self.directed, 0, 3, self.seed)
        self.assertEqual(twoPossibleWays, self.directedTwoWays)

    def testGraph(self):
        oneWay = randomWalk(self.g, 1,2)
        self.assertEqual(oneWay, self.oneWay)
        empty = randomWalk(self.g, 2, 0)
        self.assertEqual(empty, [])
        twoPossibleWays = randomWalk(self.g, 5, 7, self.seed)
        self.assertEqual(twoPossibleWays, self.twoWays)

class VertexWalkToEdgesWalkTest(unittest.TestCase):
    def setUp(self):
        self.directed = Graph([(0,1), (0,6), (1,4), (1,3), (2,3), (3,1), (1,2), (5,1)], directed=True)
        self.directedVertexPath = [0, 1, 2, 3, 1, 4]
        self.directedEdgePath = [0, 6, 4, 5, 2]
        self.directedEmptyEdge = [1]

        self.g = Graph([(0,1), (1,2), (2,3), (2,5), (3,4), (4,5)])
        self.VertexPath  = [0, 1, 2, 5, 4, 5, 2]
        self.edgePath = [0, 1, 3, 5, 5, 3]
        self.emptyPath = [5]

    def testDirectedGraph(self):
        edgePath = vertexWalkToEdgesWalk(self.directed, self.directedVertexPath)
        self.assertEqual(edgePath, self.directedEdgePath)
        empty = vertexWalkToEdgesWalk(self.directed, self.directedEmptyEdge)
        self.assertEqual(empty, [])
    
    def testGraph(self):
        edgePath = vertexWalkToEdgesWalk(self.g, self.VertexPath)
        self.assertEqual(edgePath, self.edgePath)
        empty = vertexWalkToEdgesWalk(self.directed, self.emptyPath)
        self.assertEqual(empty, [])    

class RandomWalkBetweennessTest(unittest.TestCase):
    def setUp(self):
        self.g = Graph([(0,1), (0,4), (1,4), (4,2)])
        self.vertices = [15, 15, 8, 0, 19]
        self.edges = [13, 11, 11, 10]
        self.seed = 1

    def testVertices(self):
        vertices = randomWalkBetweenness(self.g, seed = self.seed)
        self.assertTrue(np.array_equal(self.vertices, vertices))
    
    def testEdges(self):
        edges = randomWalkBetweenness(self.g, edge=True, seed = self.seed)
        self.assertTrue(np.array_equal(self.edges, edges))

class CriticalityTest(unittest.TestCase):
    def setUp(self):
        self.seed = 1
        self.delta = 0.01
        self.g = Graph([(0,1), (0,4), (1,4), (4,2)])
        self.vertexWeights = [0, 3, 4, 5, 6]
        self.edgesWeights = [2, 3, 0, 5]
        self.vertexWeightZero = 0
        self.vertexBetweennessZero = 3
        self.vertex = 1
        self.edgeWeightZero = 2
        self.edge = 3


    def testVertices(self):            
        vertexBetweenness = randomWalkBetweenness(self.g, seed=self.seed)
        #Set weights
        self.g.vs['weight'] = self.vertexWeights

        #Vertex weight 0
        expected = 0
        actual = criticality(self.g, v=self.vertexWeightZero, w='weight', s=self.seed)
        self.assertEqual(expected, actual)

        #randomWalkBerweenness 0
        expected = 0
        actual = criticality(self.g, v=self.vertexBetweennessZero, w='weight', s=self.seed)
        self.assertEqual(expected, actual)

        #General case
        expected = vertexBetweenness[self.vertex] / self.g.vs[self.vertex].attributes()['weight']
        actual = criticality(self.g, v=self.vertex, w='weight', s=self.seed)
        self.assertTrue(actual - self.delta <= expected)
        self.assertTrue(actual + self.delta >= expected)
    
    def testEdges(self):
        edgeBetweenness = randomWalkBetweenness(self.g, edge=True, seed=self.seed)
        #Set weights
        self.g.es['weight'] = self.edgesWeights

        #Edge weight 0
        expected = 0
        actual = criticality(self.g, v=self.edgeWeightZero, w='weight', edge=True, s=self.seed)
        self.assertEqual(expected, actual)

        #General case
        expected = edgeBetweenness[self.edge] / self.g.es[self.edge].attributes()['weight']
        actual = criticality(self.g, v=self.edge, w='weight', edge=True, s=self.seed)
        self.assertTrue(actual - self.delta <= expected)
        self.assertTrue(actual + self.delta >= expected)

class EntropyRankFromMatrixTest(unittest.TestCase):
    def setUp(self):
        self.m = np.array([[0,1,1], [1,0,1], [0,1,0]])
        self.vertex = 1
        self.delta = 0.1
        self.result = 0.59
    
    def testMatrix(self):
        result= entropyRankFromMatrix(self.m, self.vertex)
        self.assertTrue(result- self.delta <= self.result)
        self.assertTrue(result + self.delta >= self.result)

class EntropyRankTest(unittest.TestCase):
    def setUp(self):
        self.g = Graph([(0,1), (0,2), (1,0), (1,2), (2,1)], directed=True)
        self.vertex = 1
        self.delta = 0.1
        self.result = 0.59
    
    def testGraph(self):
        result= entropyRank(self.g, i=self.vertex)
        self.assertTrue(result- self.delta <= self.result)
        self.assertTrue(result + self.delta >= self.result)

class FreeEnergyRankTest(unittest.TestCase):
    def setUp(self):
        self.g = Graph([(0,1), (0,2), (1,0), (1,2), (2,1)], directed=True)
        self.e = 0.01
        self.vertex = 1
        self.delta = 0.1
        self.result = 0.59
    
    def testGraph(self):
        result = freeEnergyRank(self.g, i = self.vertex, e= self.e)
        self.assertTrue(result- self.delta <= self.result)
        self.assertTrue(result + self.delta >= self.result)

class BridgenessTest(unittest.TestCase):
    def setUp(self):
        self.g = Graph([(0,1), (0,3), (1,3), (1,2), (2,4), (2,5), (2,6), (4,6),(4,5), (5,6)])
        self.directed = Graph([(0,1), (0,3), (1,3), (1,2), (2,4), (2,5), (2,6), (4,6),(4,5), (5,6)], directed=True)

    def testGraph(self):
        v1 = bridgeness(self.g, i=1, j=2)
        self.assertTrue(v1 > 1.732 and v1 < 1.733, "Error in graph, case existing edge")
        self.assertRaises(Exception, bridgeness, (self.g, 0, 4))

    def testDirectedGraph(self):
        v1 = bridgeness(self.directed, i=1, j=2)
        self.assertTrue(v1 > 1.732 and v1 < 1.733, "Error in directed graph, case existing edge")
        self.assertRaises(Exception, bridgeness, (self.directed, 0, 4))

class mcvTest(unittest.TestCase):
    def setUp(self):
        self.g = Graph([(0,1), (1,2), (1,3), (1,5), (2,4), (3,4)])
        self.m = np.array(self.g.get_adjacency().data)
        self.gResult = [[0,2,3,5], [1,4], [0,4,5]]

        self.emptyGraph = Graph()
        self.emptyM = np.array(self.emptyGraph.get_adjacency().data)

        self.directed = Graph([(0,1), (1,2), (4,2), (1,5), (3,1), (3,4)], directed=True)
        self.directedMatrix = np.array(self.directed.get_adjacency().data)
        self.directedResult = [[1,2,4,5]]
        

    def testGraph(self):
        result = mcv(self.m, [], [], self.g.is_directed())
        self.assertEqual(len(result), len(self.gResult), "Error in graph, wrong number of minimal vertex covers")
        for cover in self.gResult:
            self.assertTrue(cover in result, "Error in graph, missing cover")

    def testEmptyGraph(self):
        result = mcv(self.emptyM, [], [], self.emptyGraph.is_directed())
        self.assertEqual(len(result), 0, "Error in empty graph")

    def testDirectedGraph(self):
        result = mcv(self.directedMatrix, [], [], directed =True)
        self.assertEqual(len(result), len(self.directedResult), "Error in directed graph, wrong number of minimal vertex covers")
        for cover in self.directedResult:
            self.assertTrue(cover in result, "Error in directed graph, cover not found")

class MCVTest(unittest.TestCase):
    def setUp(self):
        #One minimun cover
        self.gMinCoverOne = Graph([(0,1),(1,2)])       
        self.gMinCoverOneVertex = Graph([(0,0)])
        self.gMinCover = Graph([(0,1), (1,2), (1,3), (1,5), (2,4), (3,4)])

        #More than one minimun cover
        self.gMinCoverTwo = Graph([(0,1)])
        self.gMinCoverTwoResult = [[0], [1]]

        #Two covers of size two
        self.gMinCoverSizeTwo = Graph([(0,1), (1,2), (2,3), (3,0)])
        self.gMinCoverSizeTwoResult = [[0,2], [1,3]]


        #Directed graphs
        #One minimun cover of size one
        self.directedOneSizeOne = Graph([(0,1), (2,1)], directed=True)
        self.directedOneSizeOneResult= [[1]]

        #One of size greater than one
        self.directedOneSizeTwo = Graph([(1,0),(1,2)], directed=True)
        self.directedOneSizeTwoResult = [[0,2]]


    def testGraph(self):
        self.assertEqual(MCV(self.gMinCoverOne), [[1]], "Error in graph, must be one cover of one vertex")
        self.assertEqual(MCV(self.gMinCoverOneVertex), [[0]], "Error in graph with one vertex, must be one cover of one vertex")
        self.assertEqual(MCV(self.gMinCover), [[1,4]], "Error in graph, must be one cover of two vertex")

        result = MCV(self.gMinCoverTwo)
        self.assertEqual(len(result), len(self.gMinCoverTwoResult), "Error wrong number of found covers")
        for cover in self.gMinCoverTwoResult:
            self.assertTrue(cover in result, "Minimun cover not found")

        result2 = MCV(self.gMinCoverSizeTwo)
        self.assertEqual(len(result2), len(self.gMinCoverSizeTwoResult), "Error wrong number of found covers")
        for cover in self.gMinCoverSizeTwoResult:
            self.assertTrue(cover in result2, "Error in graph, cover not found")

    def testDirected(self):
        self.assertEqual(MCV(self.directedOneSizeOne), self.directedOneSizeOneResult, "Error in directed graph, one cover of size one")
        self.assertEqual(MCV(self.directedOneSizeTwo), self.directedOneSizeTwoResult, "Error in directed graph, one cover of size two")


class coveringDegreeTest(unittest.TestCase):
    def setUp(self):
        self.g = Graph([(0,1), (1,2), (1,5), (2,3),(5,3)])
        self.gResult0 = 2
        self.gResult1 = 1
        self.gResult4 = 0
        self.directed = Graph([(0,1), (2,1), (0,4), (2,5), (0,3),(5,3),(5,4)], directed=True)
        self.dResult2 = 0
        self.dResult1 = 1

    def testGraph(self):
        self.assertEqual(coveringDegree(self.g, v=0), self.gResult0)
        self.assertEqual(coveringDegree(self.g, v=1), self.gResult1)
        self.assertEqual(coveringDegree(self.g, v=4), self.gResult4)
    
    def testDirectedGraph(self):
        self.assertEqual(coveringDegree(self.directed, v=2), self.dResult2)
        self.assertEqual(coveringDegree(self.directed, v=1), self.dResult1)

class coveringIndexTest(unittest.TestCase):
    def setUp(self):
        self.g = Graph([(0,1), (1,2), (1,5), (2,3),(5,3)])
        self.gResult4 = 0
        self.gResult2 = 1/3
        self.e = 0.001
        self.directed = Graph([(0,1), (2,1), (0,4), (2,5), (0,3),(5,3),(5,4)], directed=True)
        self.dResult1 = 2
        self.dResult2 = 0

    
    def testGraph(self):
        self.assertEqual(coveringIndex(self.g, v=4), self.gResult4)
        self.assertTrue(coveringIndex(self.g, v=2) > self.gResult2 - self.e)
        self.assertTrue(coveringIndex(self.g, v=2) < self.gResult2 + self.e)

    def testDirectedGraph(self):
        self.assertEqual(coveringIndex(self.directed, v=1), self.dResult1)
        self.assertEqual(coveringIndex(self.directed, v=2), self.dResult2)

class weightedMatrixTest(unittest.TestCase):
    def setUp(self):
        self.g = Graph([(0,6), (0,5), (0,1), (1,2), (2,3), (1,4), (4,5), (1,5)])
        self.g.es['weight'] = [7, 6, 1, 2, 10, 2, 3, 4]
        self.weighted = np.array([[0,1,0,0,0,6,7], [1,0,2,0,2,4,0], [0,2,0,10,0,0,0], [0,0,10,0,0,0,0], [0,2,0,0,0,3,0], [6,4,0,0,3,0,0], [7,0,0,0,0,0,0]])

    def test(self):
        m = weightedMatrix(self.g, 'weight')
        self.assertTrue(np.array_equal(self.weighted, m))

if __name__ == '__main__':
    unittest.main()