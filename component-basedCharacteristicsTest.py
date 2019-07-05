from igraph import *
import unittest
from component-basedCharacteristics import *

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

if __name__ == '__main__':
    unittest.main()