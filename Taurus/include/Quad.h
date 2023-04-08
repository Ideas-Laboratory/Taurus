#ifndef _QUAD_H_
#define _QUAD_H_
#include "Node.h"
#include "iostream"
#include "math.h"
#include "vector"

#ifndef NULL
#define NULL 0
#endif
    class Quad {
    public:
        bool is_leaf;
        float size;
        float mass[2];
        float w;
        float h;
        bool have_node;
        int nodeindex;
        float center[2];
        Quad **children;

        void add(float x, float y, int index);

        //void print();
        //void addDebug(float x, float y, int index);
        Quad(std::vector<Node> &nodes, int N);

        Quad(float cx, float cy, float width, float height);

        ~Quad();
    };
    void QuadForce(Quad* qt, float& fx, float& fy, float posx, float posy, float theta, float alpha, int nodeindex,float dp,int layer=0);   
    void Quad_t_Force(Quad* qt, float& fx, float& fy, float posx, float posy, float theta, float alpha, int nodeindex,float dp,float tf,int layer=0);  
#endif
