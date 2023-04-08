## Expand Taurus with FFT solver

#### 	The current version is consistent with the features and results in the VIS22 paper.  It means code doesn't include the FFT solver.

#### How to use Taurus
```shell
cd build
./compile.sh
./main filename
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
