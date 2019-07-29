# robustness-and-resilience-in-graphs

Requirements:
- python 3.6
- igraph 0.7.1
- numpy
- scipy
- sympy

4.1 Characteristics for ranking elements

- [x] Vertex load + test
- [x] Betweenness (Already in igraph)
- [x] Random-walk betweenness + test
- [x] Criticality + test
- [x] Entropy rank + test
- [x] Free energy rank + test
- [x] Bridgeness + test
- [x] Covering degree + test
- [x] Covering index + test
- [x] Sensitivity

Auxiliary features
- [x] Random walk + test
- [x] Vertex walk to edges walk + test
- [x] Entropy rank from matrix + test
- [x] mcv(G): The set of minimal vertex covers of G + test
- [x] MCV(G): The set of minimun vertex covers of G + test
- [x] weightedMatrix(g, w): Weighted matrix of graph g, where w represents the weight attribute name in the graph + test

5.1 Component-based characteristics

- [x] Splitting number
- [x] Random-robustness index
- [x] Robustness measure + test
- [x] Connectivity robustness function + test
- [x] k-resilence factor + test
- [x] Resilence factor + test
- [x] Perturbation score
- [x] Preferential perturbation
- [x] Maximun perturbation score

Auxiliary features:
- [x] sizeMaxComponent
- [x] perturbationScoreTwo

5.2 Path-based characteristic
- [x] Pairwise disconnectivity index + test
- [x] Fragmentation
- [x] Self-sufficiency + test
- [x] k vertex-failure resilience + test
- [x] Vertex resilience + test
- [x] k edge-failure resilience + test
- [x] Egde resilience + test
- [x] Path diversity + test
- [x] Percolated path
- [x] Percolation centrality

Auxiliary features
- [x] getSimplePath(G, s, d) + test
- [x] resilience

5.3 Degree-based characteristics

- [x] Degree entropy
- [x] Relative entropy

Auxiliary features
- [x] getProbabilityDegree : The probability p(k) of the degree distribution
- [x] getDegreeDistribution

5.4 Density-based characteristics

- [x] Hub density
- [x] No name, definition523

Auxiliary features
- [x] Max degree product
- [x] Max edge betweenness

5.5 Distance-based characteristics

- [x] Geographical diversity
- [x] Effective geographical path diversity
- [x] Total graph geographical diversity
- [x] Compensated total graph geographical diversity
- [x] Functionality
- [x] Functionality loss
- [x] Global functionality loss
- [x] Temporal efficiency
- [x] Delta efficiency
- [x] Fragility
- [x] Dynamic robustness

Auxiliary features
- [x] Shortest temporal distance
- [x] Get all simple paths between two vertices
- [x] Get list of edges connected to vertex k

5.6 Random-walk based characteristics

- [x] Network criticality

Auxiliary features
- [x] Criticality of vertex
- [x] Criticality of edge

5.7 Flow-based characteristics

- [x] Electrical nodal robustness
- [x] Relative area index

5.8 Spectral characteristics

- [x] Reconstructability coefficient
- [x] Normalized subgraph centrality
- [x] Generalized robustness index
- [x] Local neutral connectivity
- [x] Closed-walk number
- [x] Redundancy of alternative paths
- [x] Neutral connectivity
- [x] Subgraph centrality
- [x] Normalized local natural connectivity

Auxiliary features
- [x] sortEigenValuesVectors
- [x] naturalConAux : Natural connectivity auxiliary function

5.9 Other characteristics

- [x] Effective graph resistence
- [x] Viral conductance

Auxiliary features
- [x] probis(g, i, s): probability that node i is infected at steady state s
- [x] y(g, s): Fraction of infected nodes