import networkx as nx
import matplotlib as mpl
mpl.use('Agg') # disable pylab's attempt to open a window
from matplotlib.pylab import figure, savefig

examples = [(nx.generators.random_graphs.barabasi_albert_graph(15, 3), 'nxba'),
            (nx.generators.classic.barbell_graph(6, 2), 'nxbarbell'),
            (nx.generators.classic.circular_ladder_graph(12), 'nxcirc'),
            (nx.generators.classic.complete_graph(12), 'nxcomplete'),
            (nx.generators.random_graphs.gnm_random_graph(20, 30), 'nxerdos'),
            (nx.generators.lattice.hexagonal_lattice_graph(3, 3), 'nxhex'),
            (nx.generators.classic.ladder_graph(12), 'nxladder'),
            (nx.generators.classic.path_graph(12), 'nxpath'),
            (nx.generators.random_graphs.random_regular_graph(3, 16), 'nxregular'),
            (nx.generators.classic.star_graph(20), 'nxstar'),
            (nx.generators.lattice.grid_2d_graph(5, 4), 'nxsquare'),
            (nx.generators.lattice.triangular_lattice_graph(5, 4), 'nxtri'),
            (nx.generators.random_graphs.connected_watts_strogatz_graph(20, 4, 0.15), 'nxws'),
            (nx.generators.random_graphs.powerlaw_cluster_graph(20, 4, 0.15), 'nxtree'),
            (nx.generators.classic.wheel_graph(12), 'nxwheel')]

redraw_examples = False
if redraw_examples:
    for (G, label) in examples:
        fig = figure()
        nx.draw(G)
        savefig(label + '.png', width = 1000)

models = [(nx.generators.random_graphs.barabasi_albert_graph, -2),
            (nx.generators.classic.barbell_graph, 2),
            (nx.generators.classic.circular_ladder_graph, 1),
            (nx.generators.classic.complete_graph, 1),
            (nx.generators.random_graphs.gnm_random_graph, 2),
            (nx.generators.lattice.hexagonal_lattice_graph, 2),
            (nx.generators.classic.ladder_graph, 1),
            (nx.generators.classic.path_graph, 1),
            (nx.generators.random_graphs.random_regular_graph, -1),
            (nx.generators.classic.star_graph, 1),
            (nx.generators.lattice.grid_2d_graph, 2),
            (nx.generators.lattice.triangular_lattice_graph, 2),
            (nx.generators.random_graphs.connected_watts_strogatz_graph, 3),
            (nx.generators.random_graphs.powerlaw_cluster_graph, 3),
            (nx.generators.classic.wheel_graph, 1)]

from math import floor, ceil, sqrt, log

def create(m, n, pr, k = None, h = None):
    G = None
    try:
        if pr == 1:
            G =  m(n)
        elif pr == 2:
            if k * k % 2 == 0:
                G = m(k, k)
            else:
                G = m(k, k + 1)                    
        elif pr == 3: # only the WS
            G =  m(n, k // 2, log(1 + (11- power) / 50))
        elif pr == -2: # BA
            G = m(n, h)
        elif pr == -1: # regular
            G =  m(h, n)                
    except Exception as e:
        print('# ERROR', m.__name__, 'failed to create a graph:', e)        
        pass
    if G is None:
        print('# OMITTING', m.__name__, 'as it produced a nil output with', n, pr, k, h)
        return None
    CC = sorted(nx.connected_components(G), key=len, reverse=True)
    if len(G) > len(CC[0]):
        print('# WARNING', m.__name__, 'resulted in a disconnected graph with', n, pr, k, h)        
        return G.subgraph(CC[0])
    else:
        return G

# <import.lst>
from characteristicsForRankingElements import vertexLoad
from characteristicsForRankingElements import randomWalkBetweenness
from characteristicsForRankingElements import criticality
from characteristicsForRankingElements import entropyRank
from characteristicsForRankingElements import freeEnergyRank
from characteristicsForRankingElements import bridgeness
from characteristicsForRankingElements import coveringDegree
from characteristicsForRankingElements import coveringIndex
from characteristicsForRankingElements import sensitivity
from componentBasedCharacteristics import splittingNumber
from componentBasedCharacteristics import randomRobustnessIndex
from componentBasedCharacteristics import robustnessMeasure53
from componentBasedCharacteristics import connectivityRobustnessFunction
from componentBasedCharacteristics import kResilienceFactor
from componentBasedCharacteristics import resilienceFactor
from componentBasedCharacteristics import perturbationScore
from degreeBasedCharacteristics import degreeEntropy
from degreeBasedCharacteristics import relativeEntropy
from densityBasedCharacteristics import hubDensity
from densityBasedCharacteristics import definition523
from distanceBasedCharacteristics import geographicalDiversity
from distanceBasedCharacteristics import effectiveGeographicalPathDiversity
from distanceBasedCharacteristics import totalGraphGeographicalDiversity
from distanceBasedCharacteristics import compensatedTotalGeographicalGraphDiversity
from distanceBasedCharacteristics import functionality
from distanceBasedCharacteristics import functionalityLoss
from distanceBasedCharacteristics import globalFunctionalityLoss
from distanceBasedCharacteristics import temporalEfficiency
from distanceBasedCharacteristics import deltaEfficiency
from distanceBasedCharacteristics import fragility
from distanceBasedCharacteristics import dynamicFragility
from flowBasedCharacteristics import electricalNodalRobustness
from flowBasedCharacteristics import relativeAreaIndex
from otherCharacteristics import effectiveGraphResistance
from otherCharacteristics import viralConductance
from pathBasedCharacteristics import pairwiseDisconnectivityIndex
from pathBasedCharacteristics import fragmentation
from pathBasedCharacteristics import selfSufficiency
from pathBasedCharacteristics import kVertexFailureResilience
from pathBasedCharacteristics import vertexResilience
from pathBasedCharacteristics import kEdgeFailureResilience
from pathBasedCharacteristics import edgeResilience
from pathBasedCharacteristics import pathDiversity
from pathBasedCharacteristics import percolatedPath
from pathBasedCharacteristics import percolationCentrality
from randomWalkBasedCharacteristics import networkCriticality
from spectralCharacteristics import reconstructabilityCoefficient
from spectralCharacteristics import normalizedSubgraphCentrality
from spectralCharacteristics import generalizedRobustnessIndex
from spectralCharacteristics import localNaturalConnectivity
from spectralCharacteristics import closedWalkNumber
from spectralCharacteristics import redundancyOfAlternativePaths
from spectralCharacteristics import naturalConnectivity
from spectralCharacteristics import subgraphCentrality
from spectralCharacteristics import normalizedLocalNaturalConnectivity
# </import.lst>

characteristics = [vertexLoad,
randomWalkBetweenness,
criticality,
entropyRank,
freeEnergyRank,
bridgeness,
coveringDegree,
coveringIndex,
sensitivity,
splittingNumber,
randomRobustnessIndex,
robustnessMeasure53,
connectivityRobustnessFunction,
kResilienceFactor,
resilienceFactor,
perturbationScore,
degreeEntropy,
relativeEntropy,
hubDensity,
definition523,
geographicalDiversity,
effectiveGeographicalPathDiversity,
totalGraphGeographicalDiversity,
compensatedTotalGeographicalGraphDiversity,
functionality,
functionalityLoss,
globalFunctionalityLoss,
temporalEfficiency,
deltaEfficiency,
fragility,
dynamicFragility,
electricalNodalRobustness,
relativeAreaIndex,
effectiveGraphResistance,
viralConductance,
pairwiseDisconnectivityIndex,
fragmentation,
selfSufficiency,
kVertexFailureResilience,
vertexResilience,
kEdgeFailureResilience,
edgeResilience,
pathDiversity,
percolatedPath,
percolationCentrality,
networkCriticality,
reconstructabilityCoefficient,
normalizedSubgraphCentrality,
generalizedRobustnessIndex,
localNaturalConnectivity,
closedWalkNumber,
redundancyOfAlternativePaths,
naturalConnectivity,
subgraphCentrality,
normalizedLocalNaturalConnectivity]
    
permitted = 2 # seconds until the execution of an individual measurement is terminated
import igraph as ig 
from time import time
from multiprocessing import Process
from math import floor, ceil, sqrt, log

def measure(c, descr, g, g2 = None):
    before, after, value = None, None, None
    if g2 == None:
        try:
            before, value, after = time(), c(g), time()
        except Exception as e:
            print('# ERROR measuring', descr, e)
            return
    else:
        try:
            before, value, after = time(), c(g, g2), time()
        except:
            print('# ERROR measuring', descr, e)
            return            
    if value is not None:
        if hasattr(value, "__iter__"):
            value = '{:f} (avg)'.format(sum(value) / len(value))
        print('{:s} {:s} {:f}'.format(descr, str(value), after - before))

for power in range(5, 10):
    n = 2**power
    k = int(floor(sqrt(n)))
    h = int(ceil(sqrt(k)))
    for (m, p) in models:
        for r in range(10 - power // 2): # less replicas for larger graphs
            G = create(m, n, p, k, h)
            if G is not None:
                G = nx.convert_node_labels_to_integers(G)
                for vertex in G:
                    if 'pos' in G.nodes[vertex]:
                        del G.nodes[vertex]['pos']
                nx.write_graphml(G, 'G.graphml')
                g = ig.read('G.graphml', format="graphml")
                d = '{:s} {:d} '.format(m.__name__, n)                
                for c in characteristics:
                    proc = Process(target=measure, args=(c, d + c.__name__, g))
                    proc.start()
                    proc.join(permitted)
                    if proc.is_alive():
                        proc.terminate()
                        proc.join()

# <import2.lst>
from componentBasedCharacteristics import preferentialPerturbation
from componentBasedCharacteristics import maximumPerturbationScore
# </import2.lst>

comparative = [preferentialPerturbation, maximumPerturbationScore]

for power in range(4, 9, 2):
    n = 2**power
    k = int(floor(sqrt(n)))
    h = int(ceil(sqrt(k)))
    for (m1, p1) in models:
        for r in range(10 - power // 2): # less replicas for larger graphs
            G1 = create(m1, n, p1, k, h)
            if G1 is not None:
                G1 = nx.convert_node_labels_to_integers(G1)
                for vertex in G1:
                    if 'pos' in G1.nodes[vertex]:
                        del G1.nodes[vertex]['pos']
                nx.write_graphml(G1, 'G1.graphml')
                g1 = ig.read('G1.graphml', format="graphml")
                for (m2, p2) in models:
                    G2 = create(m2, n, p2, k, h)
                    if G2 is not None:
                        G2 = nx.convert_node_labels_to_integers(G2)
                        for vertex in G2:
                            if 'pos' in G2.nodes[vertex]:
                                del G2.nodes[vertex]['pos']
                        nx.write_graphml(G2, 'G2.graphml')
                        g2 = ig.read('G2.graphml', format="graphml")
                        d = '{:s} {:s} {:d} '.format(m1.__name__ , m2.__name__, n)
                        for c in comparative:
                            p = Process(target=measure, args=(c, d + c.__name__, g1, g2))
                            p.start()
                            p.join(permitted)
                            if p.is_alive():
                                p.terminate()
                                p.join()
                                
                            
