<!-- # taurus.github.io -->
<!-- # taurusgl.github.io -->
<!-- # taurusgl.github.io -->
[Taurus](https://ideas-laboratory.github.io/Taurus/): Towards a Unified Force Representation and Universal Solver for Graph Layout
====
Taurus is a general graph layout framework to unify popular graph layout methods. The novelty of Taurus consists of two major components: a unified quotient-based force representation to model repulsive and attractive forces of different graph drawing techniques, and a universal augmented stochastic gradient descent (SGD) solver to find the optimal graph layout results. We release this C++ graph layout library based on Taurus that facilitates convenient implementation of graph visualizations in a unified manner as well as customizing one's own graph layout techniques.

### Example
Taurus ***online demo*** for this library can be found [here](https://ideas-laboratory.github.io/Taurus/)  (or click title).  
We compiled the C++ library into a WebAssembly module and load it into the JavaScript application, thus enabling the library to run on the browser side. WebAssembly is a low-level assembly-like language that can run C/C++ code on the Web at near-native speed.

### Install
windows and linux environment  
gcc(8.10.0) cmake(3.19.2)

### Build
We have compiled Taurus into a static library. Here is a test file that demonstrates how to use the Taurus library. To run this test-file, run build.sh and the file will be compiled and linked to the lib.
```bash
./build.sh
./test
```

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
Taurus implements six graph layout methods:***SM, FDP, LinLog, Maxent, FMMM, BSM***. You can use these methods by inputting a graph by calling the corresponding interface.
```bash
SMLayout(g)
FDPLayout(g)
...
```
The library also supports customizing one's own graph layout. You can set the force model and pass it to ***Layout*** interface.
```bash
Layout(g,force)
```
Layout results are stored in *g* as coordinates. The library provides an interface to draw the layout results as svg.
```bash
g.drawSVG(filename,method-name)
```
### Licensing
The source code is licensed under GPL v3. License is available [here](/LICENSE).

If you have problems, please submit an issue.
