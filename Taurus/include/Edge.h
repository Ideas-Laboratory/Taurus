#ifndef _EDGE_H_
#define _EDGE_H_
 typedef struct Edge {
        int w;
        int end;
        int edge_strength;
        float edge_redundancy;

        bool operator<(const Edge &e) const {
            return edge_strength < e.edge_strength;
        }
    } Edge;
#endif