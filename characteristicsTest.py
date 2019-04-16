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
        self.mcvg = [[0,2,3,5], [1,4], [0,4,5]]

        self.emptyGraph = Graph()
        self.emptyM = np.array(self.emptyGraph.get_adjacency().data)

        self.directed = Graph([(0,1), (1,2), (4,2), (1,5), (3,1), (3,4)], directed=True)
        

    def testGraph(self):
        result = mcv(self.m, [], [], self.g.is_directed())
        self.assertEqual(len(result), len(self.mcvg), "Error in graph, wrong number of minimal vertex covers")
        for cover in self.mcvg:
            self.assertTrue(cover in result, "Error in graph, missing cover")

    def testEmptyGraph(self):
        result = mcv(self.emptyM, [], [], self.emptyGraph.is_directed())
        self.assertEqual(len(result), 0, "Error in empty graph")

    def testDirectedGraph(self):



class CoveringDegreeTest(unittest.TestCase):
    def setUp(self):
        self.directed = Graph([(0,1), (0,2), (1,2), (1,3), (2,3)], directed=True)
        self.g = Graph([(0,1), (1,2), (1,3), (1,5), (2,4), (3,4)])

    def testDirectedGraph(self):
        self.assertEqual(coveringDegree(self.directed, 0), 0, "Error in directed graph, case v has no incident edges")
        self.assertEqual(coveringDegree(self.directed, 1), 1, "Error in directed graph, case v has incident edges")
        #TODO: Buscar ejemplo para grafo dirigido con resultado > 1

    def testGraph(self):
        self.assertEqual(coveringDegree(self.g, 5), 2)

    
if __name__ == '__main__':
    unittest.main()