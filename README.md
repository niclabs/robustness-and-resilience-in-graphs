# robustness-and-resilience-in-graphs

Requirements:
- python3
- igraph https://igraph.org/c/doc/igraph-Flows.html
- numpy
- scipy
- sympy
- networkx

Executing the experiments additionally requires networkx, and the generation of the example figures employs matplotlib.

4.1 Characteristics for ranking elements

- [x] Vertex load
- [x] Betweenness (Already in igraph)
- [x] Random-walk betweenness
- [x] Criticality, graph needs a 'weight' attribute on nodes or edges as appropriate
- [x] Entropy rank
- [x] Free energy rank
- [x] Bridgeness
- [x] Covering degree
- [x] Covering index
- [x] Sensitivity, graph needs a 'weight' attribute on edges

5.1 Component-based characteristics
- [x] Normalized GCC
- [x] Splitting number
- [x] Random-robustness index
- [x] Robustness measure
- [x] Connectivity robustness function
- [x] k-resilience factor
- [x] Resilience factor
- [x] Perturbation score
- [x] Preferential perturbation
- [x] Maximun perturbation score
- [x] Robustness index
- [x] Robustness measure R 

5.2 Path-based characteristic
- [x] Local connectivity
- [x] Global connectivity
- [x] Pairwise disconnectivity index
- [x] Fragmentation
- [x] Self-sufficiency
- [x] k vertex-failure resilience
- [x] Vertex resilience
- [x] k edge-failure resilience
- [x] Egde resilience
- [x] Path diversity
- [x] Percolated path, graph needs attribute 'state' on nodes
- [x] Percolation centrality, graph needs attribute 'state' on nodes
- [x] Percentage of noncritical nodes
- [x] Treeness

5.3 Degree-based characteristics

- [x] Degree entropy
- [x] Relative entropy

5.4 Density-based characteristics

- [x] Hub density
- [x] No name, definition523

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
- [x] Vulnerability

5.6 Random-walk based characteristics

- [x] Network criticality, graph needs attribute 'weight' on edges

5.7 Flow-based characteristics

- [x] Electrical nodal robustness, graphs needs a 'flow' attribute on edges
- [x] Relative area index, return None when a value is not valid

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

5.9 Other characteristics
- [x] Entropy
- [x] Effective graph resistence, graph needs 'weight' attribute on edges
- [x] Viral conductance
- [x] RCB