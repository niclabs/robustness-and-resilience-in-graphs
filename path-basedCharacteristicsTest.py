from igraph import *
import unittest
from path-basedCharacteristics import *

class pairwiseDisconnectivityIndexTest(unittest.TestCase):
    def setUp(self):
        self.graph = Graph([(0,1), (0,2), (1,5), (2,3), (3,4), (4,5)])
        self.graphv2 = 2
        self.graphv2Result = 1/3
        self.delta = 0.01

        self.directed = Graph([(0,2), (0,4), (1,5), (2,1), (3,1), (5,3)], directed=True)
        self.directedv1 = 1
        self.directedv1Result = 11/14
    
    def testGraph(self):
        r1 = pairwiseDisconnectivityIndex(self.graph, self.graphv2)
        self.assertTrue(r1 > self.graphv2Result - self.delta)
        self.assertTrue(r1 < self.graphv2Result + self.delta)
    
    def testDirectedGraph(self):
        r1 = pairwiseDisconnectivityIndex(self.directed, self.directedv1)
        self.assertTrue(r1 > self.directedv1Result - self.delta)
        self.assertTrue(r1 < self.directedv1Result + self.delta)

class selfSufficiencyTest(unittest.TestCase):
    def setUp(self):
        self.simpleGraph = Graph([(0,3), (1,2), (2,3)])
        self.ssListSimple = [[[0, 1, 2], [3]], [[1, 3], [2]], [[2, 0], [1]], [[1, 3, 2], [0]]]
        self.nsListSimple = [[[0, 1, 2], [3]], [[1, 3], [2]], [[2, 0], [7]], [[1, 3, 2], [0]]]
        self.GraphTwoClusters = Graph([(0,1), (0,2), (0,3), (1,3), (1,2), (2,3), (4,5), (4,6), (4,7), (5,7), (5,6), (6,7)])
        self.ssListTwoClusters = [[[0, 1, 2], [3]], [[0, 1, 3], [2]], [[0, 2, 3], [1]], [[1, 2, 3], [0]], [[4, 5, 6], [7]], [[4, 5, 7], [6]], [[4, 6, 7], [5]], [[5, 6, 7], [4]]]
        self.nssListTwoClusters = [[[0, 1, 2], [3]], [[0, 1, 3], [2]], [[0, 2, 3], [1]], [[1, 2, 3], [5]], [[4, 5, 6], [7]], [[4, 5, 7], [6]], [[4, 6, 7], [5]], [[5, 6, 7], [4]]]

    def testTrueSimple(self):
        self.assertTrue(selfSufficiency(self.simpleGraph, self.ssListSimple))

    def testFalseSimple(self):
        self.assertFalse(selfSufficiency(self.simpleGraph, self.nsListSimple))

    def testTrueTwoClusters(self):
        self.assertTrue(selfSufficiency(self.GraphTwoClusters, self.ssListTwoClusters))
    
    def testFalseTwoClusters(self):
        self.assertFalse(selfSufficiency(self.GraphTwoClusters, self.nssListTwoClusters))

class kVertexFailureResilienceTest(unittest.TestCase):
    def setUp(self):
        self.simpleGraph = Graph([(0,3), (1,2), (2,3), (3,4)])
        self.ssListSimple = [[[0, 1, 2], [3]], [[1, 3], [2]], [[2, 0], [1]], [[1, 3, 2], [0]], [[3], [1]]]

        self.graph = Graph([(0,1), (0,2), (0,3), (1,3), (1,2), (2,3)])
        self.graphList = [ [[0, 1, 2], [3]], [[0, 1, 3], [2]], [[0, 2, 3], [1]], [[1, 2, 3], [0]]]
    
    def testSimpleGraph(self):
        #k = 0
        self.assertTrue(kVertexFailureResilience(self.simpleGraph, self.ssListSimple, 0))

        #k = 1
        self.assertFalse(kEdgeFailureResilience(self.simpleGraph, self.ssListSimple, 1))

        #k = 2
        self.assertFalse(kEdgeFailureResilience(self.simpleGraph, self.ssListSimple, 2))
    
    def testSimpleGraphException(self):
        with self.assertRaises(Exception) as context:
            kVertexFailureResilience(self.simpleGraph, self.ssListSimple, self.simpleGraph.vcount() + 1)

        self.assertTrue('Number of vertices to fail can not be greater than the total vertices' in str(context.exception))

    def testGraph(self):
        self.assertTrue(kVertexFailureResilience(self.graph, self.graphList, 1))
        self.assertTrue(kVertexFailureResilience(self.graph, self.graphList, 2))
        self.assertFalse(kVertexFailureResilience(self.graph, self.graphList, 3))

class vertexResilienceTest(unittest.TestCase):
    def setUp(self):
        self.graph = Graph([(0,1), (0,2), (0,3), (1,3), (1,2), (2,3)])
        self.graphList = [ [[0, 1, 2], [3]], [[0, 1, 3], [2]], [[0, 2, 3], [1]], [[1, 2, 3], [0]]]
        self.graphResult = 2

        self.nsGraph = Graph([(0,2), (0,3), (0,4), (2,4), (3,4)])
        self.nsGraphList = [[[0, 1, 2], [3]], [[1], [2]], [[3, 1], [2]], [[2, 3], [1]], [[4], [0]]]
        self.nsGraphResult = 0

    def testGraph(self):
        self.assertEqual(vertexResilience(self.graph, self.graphList), self.graphResult)

    def testnsGraph(self):
        self.assertEqual(vertexResilience(self.nsGraph, self.nsGraphList), self.nsGraphResult)

class kEdgeFailureResilienceTest(unittest.TestCase):
    def setUp(self):
        self.simpleGraph = Graph([(0,3), (1,2), (2,3), (3,4)])
        self.ssListSimple = [[[0, 1, 2], [3]], [[1, 3], [2]], [[2, 0], [1]], [[1, 3, 2], [0]], [[3], [1]]]

        self.graph = Graph([(0,1), (0,2), (0,3), (1,3), (1,2), (2,3)])
        self.graphList = [ [[0, 1, 2], [3]], [[0, 1, 3], [2]], [[0, 2, 3], [1]], [[1, 2, 3], [0]]]
    
    def testSimpleGraph(self):
        #k = 0
        self.assertTrue(kEdgeFailureResilience(self.simpleGraph, self.ssListSimple, 0))

        #k = 1
        self.assertFalse(kEdgeFailureResilience(self.simpleGraph, self.ssListSimple, 1))

    def testSimpleGraphException(self):
        with self.assertRaises(Exception) as context:
            kEdgeFailureResilience(self.simpleGraph, self.ssListSimple, self.simpleGraph.ecount() + 1)

        self.assertTrue('Number of edges to fail can not be greater than the total edges' in str(context.exception))

    def testGraph(self):
        self.assertTrue(kEdgeFailureResilience(self.graph, self.graphList, 0))
        self.assertTrue(kEdgeFailureResilience(self.graph, self.graphList, 1))
        self.assertTrue(kEdgeFailureResilience(self.graph, self.graphList, 2))
        self.assertFalse(kEdgeFailureResilience(self.graph, self.graphList, 3))

class edgeResilienceTest(unittest.TestCase):
    def setUp(self):
        self.graph = Graph([(0,1), (0,2), (0,3), (1,3), (1,2), (2,3)])
        self.graphList = [ [[0, 1, 2], [3]], [[0, 1, 3], [2]], [[0, 2, 3], [1]], [[1, 2, 3], [0]]]
        self.graphResult = 2

        self.nsGraph = Graph([(0,2), (0,3), (0,4), (2,4), (3,4)])
        self.nsGraphList = [[[0, 1, 2], [3]], [[1], [2]], [[3, 1], [2]], [[2, 3], [1]], [[4], [0]]]
        self.nsGraphResult = 0
    
    def testGraph(self):
        self.assertEqual(edgeResilience(self.graph, self.graphList), self.graphResult)

    def testnsGraph(self):
        self.assertEqual(edgeResilience(self.nsGraph, self.nsGraphList), self.nsGraphResult)


class getSimplePathTest(unittest.TestCase):
    def setUp(self):
        self.seed = 5
        self.directed = Graph([(0,1), (0,2), (1,3), (2,3), (2, 4)], directed = True)
        self.directedOneWay = [[0, 2, 4], [1, 4]]
        self.directedTwoWays = [[0, 1 , 3], [0, 2]]

        self.graph = Graph([(0,1), (1,2), (1,3), (1,5), (3,5)])
        self.graphOneWay = [[0, 1, 2], [0, 1]]
        self.graphTwoWays = [[0, 1, 3], [0, 2]]
    
    def testDirectedGraph(self):
        #There is only one way, with or without seed
        self.assertEqual(getSimplePath(self.directed, 0, 4, 0), self.directedOneWay)
        self.assertEqual(getSimplePath(self.directed, 0, 4, self.seed), self.directedOneWay)

        #There is two ways
        self.assertEqual(getSimplePath(self.directed, 0, 3, self.seed), self.directedTwoWays)

    def testDirectedGraphException(self):
        with self.assertRaises(Exception) as context:
            getSimplePath(self.directed, 3, 4, 0)

        self.assertTrue('There is no path between s and d' in str(context.exception))
    
    def testGraph(self):
        #There is only one way, with or without seed
        self.assertEqual(getSimplePath(self.graph, 0, 2, 0), self.graphOneWay)
        self.assertEqual(getSimplePath(self.graph, 0, 2, self.seed), self.graphOneWay)

        #There is two ways
        self.assertEqual(getSimplePath(self.graph, 0, 3, self.seed), self.graphTwoWays)
    
    def testgraphException(self):
        with self.assertRaises(Exception) as context:
            getSimplePath(self.graph, 2, 4, 0)

        self.assertTrue('There is no path between s and d' in str(context.exception))

class pathDiversityTest(unittest.TestCase):
    def setUp(self):
        self.graph = Graph([(0,1), (1,2), (1,3), (1,5), (3,5)])
        self.seed = 6
        self.result = 1 - (4/5)
        self.delta = 0.01

        self.directed = Graph([(0,1), (0,2), (1,3), (2,4), (2,5), (3,4)], directed = True)
        self.directedResult = 1 - (2/5)
    
    def testGraph(self):
        result= pathDiversity(self.graph, 0, 3, self.seed)
        self.assertTrue(result > self.result - self.delta)
        self.assertTrue(result < self.result + self.delta)

        self.assertEqual(pathDiversity(self.graph, 0, 2, self.seed), 0)
    
    def testDirectedGraph(self):
        result = pathDiversity(self.directed, 0, 4, self.seed)
        self.assertTrue(result > self.directedResult - self.delta)
        self.assertTrue(result < self.directedResult + self.delta)

        self.assertEqual(pathDiversity(self.graph, 0, 5, self.seed), 0)
    
    def testDirectedGraphException(self):
        with self.assertRaises(Exception) as context:
            pathDiversity(self.directed, 3, 0, 0)
        
        self.assertTrue('There is no path between s and d' in str(context.exception))
    
if __name__ == '__main__':
    unittest.main()