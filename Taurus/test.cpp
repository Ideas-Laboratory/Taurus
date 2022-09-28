#include "layout.h"
#include "Multilevel.h"
using namespace std;

int main(){
    //1.declare a graph
    graph g;
    //2.initialize a graph by edges

    g.initgraph("data/bus1138.mtx");
    // define the initial coordinates of the node
    g.initRandomPosition();// pmds  input position
    //solve for the length of the shortest path between any two pairs of nodes
    g.solveDij();
    //3.use graph layout methods or custom force representation model
    /*vector<force> f;
    f.push_back(force(0,1,1,2));
    f.push_back(force(0,-1,0,1));
    //4.solve the layout
    Layout(g,f);*/
    clock_t start=clock();
    SMLayout(g);
    clock_t end=clock();
    float timeuse = 1000 * (end - start) / CLOCKS_PER_SEC;
    printf("use time:%f ms\n", timeuse);
    //5.draw the layout
    g.drawSVG("bus1138.mtx","SM");


}