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

class robustnessMeasure53Test(unittest.TestCase):
    def setUp(self):
        self.graph = Graph([(0,2), (0,3), (0,4), (1,3), (1,4), (1,5), (2,3), (2,4)])
        self.result = 11/6
        self.delta = 0.001

    def testGraph(self):
        beforeNumberVertices = self.graph.vcount()
        self.assertTrue(robustnessMeasure53(self.graph) > self.result - self.delta)
        afterNumberVerices = self.graph.vcount()
        self.assertEqual(beforeNumberVertices, afterNumberVerices) #Check that no vertices have been removed
        self.assertTrue(robustnessMeasure53(self.graph) < (self.result + self.delta))

class connectivityRobustnessFunctionTest(unittest.TestCase):
    def setUp(self):
        self.twoVertexGraph = Graph([(0,1)])
        self.twoVertexGraphK = 1
        self.twoVertexGraphResult = 1
        self.graph = Graph([(0,1), (1,2)])
        self.k1 = 1
        self.k1result1 = 1
        self.k1result2 = 1/2
        self.k2 = 2
        self.k2result = 1
        self.delta = 0.001
    
    def testTwoVertexGraph(self):
        beforeNumberVertices = self.twoVertexGraph.vcount()
        self.assertEqual(connectivityRobustnessFunction(self.twoVertexGraph, self.twoVertexGraphK), self.twoVertexGraphResult) #Aqui
        afterNumberVertices = self.twoVertexGraph.vcount()
        self.assertEqual(beforeNumberVertices, afterNumberVertices)
    
    def testGraph(self):
        #k = 2
        self.assertEqual(connectivityRobustnessFunction(self.graph, self.k2), self.k2result)

        #k1
        result = connectivityRobustnessFunction(self.graph, self.k1)
        self.assertTrue(result == self.k1result1 or (result < self.k1result2 + self.delta and result > self.k1result2 - self.delta))

class kResilienceFactorTest(unittest.TestCase):
    def setUp(self):
        self.twoVertexGraph = Graph([(0,1)])
        self.twoVertexGraphK1 = 1
        self.twoVertexGraphResultK1 = 100
        self.twoVertexGraphK2 = 2
        self.twoVertexGraphResultK2 = 0
        self.twoVertexGraphK3 = 3
        self.twoVertexGraphResultK3 = 0

        self.graph = Graph([(1,3)])
        self.graphK2 = 2
        self.graphK2result = (2/3) * 100
        self.graphK3 = 3
        self.graphK3result1 = (1/3) * 100
        self.graphK3result2 = (2/3) * 100
        self.delta = 0.01
    
    def testTwoVertexGraph(self):
        #Remove 0 
        self.assertEqual(kResilienceFactor(self.twoVertexGraph, self.twoVertexGraphK1), self.twoVertexGraphResultK1)        
        beforeNumberVertices = self.twoVertexGraph.vcount()
        #Remove 1 vertex
        self.assertEqual(kResilienceFactor(self.twoVertexGraph, self.twoVertexGraphK2), self.twoVertexGraphResultK2)
        afterNumberVertices = self.twoVertexGraph.vcount()
        self.assertEqual(beforeNumberVertices, afterNumberVertices) 
        #Remove all vertices
        self.assertEqual(kResilienceFactor(self.twoVertexGraph, self.twoVertexGraphK3), self.twoVertexGraphResultK3)

    def testGraph(self):
        r1 = kResilienceFactor(self.graph, self.graphK2)
        self.assertTrue(r1 > self.graphK2result - self.delta)
        self.assertTrue(r1 < self.graphK2result + self.delta)

        r2 = kResilienceFactor(self.graph, self.graphK3)
        self.assertTrue(r2 > self.graphK3result1 - self.delta or r2 > self.graphK3result2 - self.delta)
        self.assertTrue(r2 < self.graphK3result1 + self.delta or r2 < self.graphK3result2 + self.delta)

class resilienceFactorTest(unittest.TestCase):
    def setUp(self):
        self.twoVertexGraph = Graph([(0,1)])
        self.twoVertexGraphResult = -1

        self.threeVertexGraph = Graph([(1,2)])
        self.threeVertexGraphResult = 50

        self.graph = Graph([(1,3)])
        self.graphR1 = np.array([(2/3)*100, (1/3)*100])
        self.graphR2 = np.array([(2/3)*100, (2/3)*100])
        self.graphResult1 = np.mean(self.graphR1)
        self.graphResult2 = np.mean(self.graphR2)
        self.delta = 0.01
    
    def testTwoVertexGraph(self):
        self.assertEqual(resilienceFactor(self.twoVertexGraph), self.twoVertexGraphResult)
    
    def testThreeVertexGraph(self):
        self.assertEqual(resilienceFactor(self.threeVertexGraph), self.threeVertexGraphResult)

    def testGraph(self):
        r = resilienceFactor(self.graph)
        self.assertTrue(r > self.graphResult1 - self.delta or r > self.graphResult2 - self.delta)
        self.assertTrue(r < self.graphResult1 + self.delta or r < self.graphResult2 + self.delta)


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