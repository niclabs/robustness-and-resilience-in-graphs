from igraph import *
import unittest
from characteristics import *

class VertexLoadTest(unittest.TestCase):
    def setUp(self):
        self.directed = Graph([(1,2), (1,3), (3,4)], directed =True)
        self.g =  Graph([(1,2), (1,3), (3,4)])
    
    def testDirectedGraph(self):
        self.assertEqual(vertexLoad(self.directed, 0, 1), 0, "Error in diredted graph, case v has no neighbors")
        self.assertEqual(vertexLoad(self.directed, 0, 0), 1, "Error in diredted graph, case v has no neighbors and n= 0")
        self.assertEqual(vertexLoad(self.directed, 0, 5), 0, "Error in diredted graph, case v has no neighbors and n= 5")
        self.assertEqual(vertexLoad(self.directed, 1, 1), 6, "Error in directed graph, case v has neighbors anmd n= 1")
        self.assertEqual(vertexLoad(self.directed, 1, 0), 1, "Error in directed graph, case v has neighbors and n= 0")
        self.assertEqual(vertexLoad(self.directed, 1, 5), 7776, "Error in directed graph, case v has neighbors and n= 6")

    def testGraph(self):
        self.assertEqual(vertexLoad(self.g, 0, 1), 0, "Error in graph, case v has no neighbors")
        self.assertEqual(vertexLoad(self.g, 0, 0), 1, "Error in graph, case v has no neighbors and n= 0")
        self.assertEqual(vertexLoad(self.g, 0, 5), 0, "Error in graph, case v has no neighbors and n= 5")
        self.assertEqual(vertexLoad(self.g, 1, 1), 6, "Error in graph, case v has neighbors anmd n= 1")
        self.assertEqual(vertexLoad(self.g, 1, 0), 1, "Error in graph, case v has neighbors and n= 0")
        self.assertEqual(vertexLoad(self.g, 1, 5), 7776, "Error in graph, case v has neighbors and n= 6")

class BridgenessTest(unittest.TestCase):
    def setUp(self):
        self.g = Graph([(0,1), (0,3), (1,3), (1,2), (2,4), (2,5), (2,6), (4,6),(4,5), (5,6)])
        self.directed = Graph([(0,1), (0,3), (1,3), (1,2), (2,4), (2,5), (2,6), (4,6),(4,5), (5,6)], directed=True)

    def testGraph(self):
        v1 = bridgeness(self.g, 1, 2)
        self.assertTrue(v1 > 1.732 and v1 < 1.733, "Error in graph, case existing edge")
        self.assertRaises(Exception, bridgeness, (self.g, 0, 4))

    def testDirectedGraph(self):
        v1 = bridgeness(self.directed, 1, 2)
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


class CoveringDegreeTest(unittest.TestCase):
    def setUp(self):
        self.directed = Graph([(0,1), (0,2), (1,2), (1,3), (2,3)], directed=True)
        self.g = Graph([(0,1), (1,2), (1,3), (1,6), (2,4), (3,4)])

    def testDirectedGraph(self):
        self.assertEqual(coveringDegree(self.directed, 0), 0, "Error in directed graph, case v has no incident edges")
        self.assertEqual(coveringDegree(self.directed, 1), 1, "Error in directed graph, case v has incident edges")

    def testGraph(self):
        self.assertEqual(coveringDegree(self.g, 5), 0, "Error in graph, v has no neighbors")
        self.assertEqual(coveringDegree(self.g, 6), 2, "Error in graph, case v has incident edges")
        self.assertEqual(coveringDegree(self.g, 0), 2, "Error in graph, case v has incident edges")

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
        self.assertEqual(coveringDegree(self.g, 0), self.gResult0)
        self.assertEqual(coveringDegree(self.g, 1), self.gResult1)
        self.assertEqual(coveringDegree(self.g, 4), self.gResult4)
    
    def testDirectedGraph(self):
        self.assertEqual(coveringDegree(self.directed, 2), self.dResult2)
        self.assertEqual(coveringDegree(self.directed, 1), self.dResult1)

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
        self.assertEqual(coveringIndex(self.g, 4), self.gResult4)
        self.assertTrue(coveringIndex(self.g, 2) > self.gResult2 - self.e)
        self.assertTrue(coveringIndex(self.g, 2) < self.gResult2 + self.e)

    def testDirectedGraph(self):
        self.assertEqual(coveringIndex(self.directed, 1), self.dResult1)
        self.assertEqual(coveringIndex(self.directed, 2), self.dResult2)

        









    
if __name__ == '__main__':
    unittest.main()