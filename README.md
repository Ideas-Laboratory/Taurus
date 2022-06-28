<!-- # taurus.github.io -->
<!-- # taurusgl.github.io -->
<!-- # taurusgl.github.io -->
[Taurus](https://taurusgl.github.io/): Towards a Unified Force Representation and Universal Solver for Graph Layout
====
Taurus is a general graph layout framework to unify popular graph layout methods. The novelty of Taurus consists of two major components: a unified quotient-based force representation to model repulsive and attractive forces of different graph drawing techniques, and a universal augmented stochastic gradient descent (SGD) solver to find the optimal graph layout results. We release this C++ graph layout library based on Taurus that facilitates convenient implementation of graph visualizations in a unified manner as well as customizing one's own graph layout techniques.

### Example
Taurus ***online demo*** for this library can be found [here](https://taurusgl.github.io/)  (or click title).  
We compiled the C++ library into a WebAssembly module and load it into the JavaScript application, thus enabling the library to run on the browser side. WebAssembly is a low-level assembly-like language that can run C/C++ code on the Web at near-native speed.

### Install
windows and linux environment  
gcc(8.10.0) cmake(3.19.2)

<!-- ### Build
 -->
### Start
To use Taurus, you just need to include this library and define a graph, and then you can construct the structure of the graph by input file or customize the topology of the graph.
```bash
graph g
g.initgraph(filename)
```

Before drawing a graph, you can choose the initial layout of the graph. The library provides two initial layouts: random layout (default)、PivotMDS layout. Besides, You can use customized coordinates from the input file as initial coordinates.
```bash
g.initRandomPosition
g.initPivotMDSPosition
g.initPosition(filename)
```
Taurus implements five graph layout methods:***SM, FDP, LinLog, Maxent, BSM***. You can use these methods by inputting a graph by calling the corresponding interface.
```bash
SMLayout(g)
FDPLayout(g)
...
```
The library also supports customizing one's own graph layout. You can set the force model and pass it to ***Layout*** interface.
```bash
Layout(g,force)
```
Layout results are stored in *g* as coordinates.  

*We will release the C++ library after paper acepted.*
