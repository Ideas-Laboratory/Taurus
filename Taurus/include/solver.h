
#ifndef TAURUS_SOLVER_H
#define TAURUS_SOLVER_H

#endif //TAURUS_SOLVER_H
#include "time.h"

#include <iomanip>
///#pragma GCC optimize(3,"Ofast","inline")
#pragma warning(disable: 4702)
#pragma warning(pop)
#include <ctime>
#include <iterator> //
#include <numeric>
#include <random>

#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstring>
#include "algorithm"

#include <queue>
#include <map>
#include <stack>

#include <time.h>
#include "sstream"
#include "set"
#include "randomkit.h"
#ifndef NULL
#define NULL 0
#endif
#define randf() (1.0f*random()/RAND_MAX)
#define Pair pair<int,int>
#define MAX 200000+100
#define INF 0x3f3f3f3f

#define DEFAULT_WEIGHT 1
#define velocityDecay  (0.8f)
typedef std::pair<int,int> node_pair;
using namespace std;

    struct Node {
        int id;
        float x, y;
        int deg;
        //Node(int i):index(i),x(0),y(0){}
    };

    struct Link {
        int target, source;

        Link() {}

        Link(int s, int t) : source(s), target(t) {}

        friend bool operator<(const Link &n1, const Link &n2) {
            if (n1.source != n2.source) {
                return n1.source < n2.source;

            } else {
                return n1.target < n2.target;
            }

        }
    };

    typedef struct Edge {
        int w;
        int end;
        int edge_strength;
        float edge_redundancy;

        bool operator<(const Edge &e) const {
            return edge_strength < e.edge_strength;
        }
    } Edge;

    struct para {
        float ka, q, a, kr, r, p, pmds;
        double w1 = 1, w2 = 1, w3 = 1;

        para() : ka(0), q(0), a(0), kr(0), r(0), p(0), pmds(0) {}

        para(float ka, float q, float a, float kr, float r, float p, float neighbor, double w1, double w2, double w3)
                : ka(ka), q(q), a(a), kr(kr), r(r), p(p), pmds(neighbor), w1(w1), w2(w2), w3(w3) {}

        para(float ka, float q, float a, float kr, float r, float p, float neighbor) : ka(ka), q(q), a(a), kr(kr), r(r),
                                                                                       p(p), pmds(neighbor) {}
    };

    struct constraint {
        int i, j;
        double d, w, wij, wji;
        para pa;

        constraint() {}

        constraint(int i, int j, double d, double w, para pa) : i(i), j(j), d(d), w(w), pa(pa) {}

        constraint(int i, int j, double d, double w) : i(i), j(j), d(d), w(w) {}

        constraint(int i, int j, double d, double wij, double wji) : i(i), j(j), d(d), wij(wij), wji(wji) {}

        constraint(int i, int j, double d) : i(i), j(j), d(d), wij(0), wji(0) {}
        /*bool operator<(const constraint &c)const{
            if(c.i!=i) return i<c.i;
            else return j<c.j;
        }*/
    };

    struct force {
        float range, k, a, b;

        force() {}

        force(float _k, float _a, float _b) : range(-1), k(_k), a(_a), b(_b) {}

        force(float _range, float _k, float _a, float _b) : range(_range), k(_k), a(_a), b(_b) {}
    };


    struct constraint_sgd {
        int i, j;
        double d, w;
        mutable vector<force> force_sgd;

        constraint_sgd() {}

        constraint_sgd(int _i, int _j) : i(_i), j(_j), d(0), w(0), force_sgd() {}

        constraint_sgd(int _i, int _j, double _d, double _w, vector<force> _force_sgd) : i(_i), j(_j), d(_d), w(_w),
                                                                                         force_sgd(_force_sgd) {}

        bool operator<(const constraint_sgd &_cs) const {
            if (i != _cs.i) {
                return i < _cs.i;
            } else {
                return j < _cs.j;
            }
        }
    };

    struct merge_cons {
        vector<force> force_m;
        bool flag = false;
        double d, w;
    };

    int str2int(string s);

    float str2float(string s);

    class graph {
    public:
        int n, m;
        vector<Node> nodes;
        vector<Link> links;

        vector<int> pivots;
        bool deg_rep = false;

        vector<vector<double>> sp;
        vector<constraint> constraints;
        vector<constraint> constraintbh;
        set<constraint_sgd> set_sgd;
        vector<constraint_sgd> constraints_sgd;
        vector<force> bh_forces;
        vector<float> count;
        vector<float> force_r;
        vector<float> schedule;
        vector<float> schedulet;
        rk_state rstate;
        para paraSGD;
        set<constraint> consBH;
        para paraBH;
        double alphaDecay = 0.034;

        graph();

        void readGraph(string filename);

        void readGraphMtx(string filename);

        void initgraph(string filename);

        void initEdge(int m);

        void initNode(int n);

        void initPosition();

        void initRandomPosition();

        void initPivotMDSPosition(int p = 100);

        void readPosition(string filename);

        void readGVPosition(string filename);

        void readLinlogPosition(string filename);

        void initPosition(string filename);

        void solveDij();


        void append_E_range(para pa);

        void append_N_range(para pa);

        void append_Np_range(para pa);

        void append_S_range(int neighbor, para pa);

        void append_range(vector<Link> nodeList, para pa);

        void append_BH_range(int neighbor, para pa);

        void append_E_range(force pa);

        void append_N_range(force pa);

        void append_Np_range(force pa);

        void append_S_range(int neighbor, force pa);

        void append_range(vector<Link> nodeList, force pa);

        void append_BH_range2(force pa);

        void preSolveSGD(double eps = 0.1, int t_max = 200, int seed = 42, float eta_max = 100.0f);

        void preSolveSGD2(double eps = 0.1, int t_max = 200, int seed = 42, float eta_max = 100.0f);

        void preSolveBH();

        void preSolveBH2();

        bool solveSGD(int iter);

        bool solveSGD2(int iter);

        bool power_law_graph();

        void solveBH(int iter);

        void solveBH2(int iter);

        void outToFile(string filename);

        void RelationOutToFile(string filename);

        void outToFile(string filename, string methodname);

        void RelationOutToFile(string filename, string methodname);

        void MtxToCsv(string filename);

        void drawSVG(string filename, string method);

    };

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
        Quad(vector<Node> &nodes, int N);

        Quad(float cx, float cy, float width, float height);

        ~Quad();
    };

    vector<vector<int> > build_graph_unweighted(int n, int m, int *I, int *J);

    void PivotMDS(vector<Node> &nodes, int k, int N, vector<vector<double>> &shortPat, vector<int> &pivot);

    void Dij(int V, int E, vector<Edge> e[MAX], int st, vector<float> &dist);

    vector<int> maxmin_random_sp_unweighted(const vector<vector<int> > &graph, int n_pivots, int p0, int seed);

    void maxmin_bfs_unweighted(const vector<vector<int> > &graph, const int p, vector<int> &mins, vector<int> &argmins);
