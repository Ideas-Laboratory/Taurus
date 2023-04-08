
#ifndef TAURUS_SOLVER_H
#define TAURUS_SOLVER_H



#include "time.h"
//#include <io.h>
//#include "direct.h"
#include <iomanip>
///#pragma GCC optimize(3,"Ofast","inline")
#pragma warning(disable: 4702) //��ֹ
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

#include "PivotMDS.h"

//#include <sys/time.h>
#include <time.h>
#include "sstream"
#include "set"
#include "shuffle.h"
//#include "randomkit.h"
//#include "Link.h"

//#include "range.h"
#include "Quad.h"

#define randf() (1.0f*random()/RAND_MAX)


#define INF 0x3f3f3f3f

#define DEFAULT_WEIGHT 1
#define velocityDecay  (0.8f)

    int str2int(std::string s);

    float str2float(std::string s);

    class graph {
    public:
        int n, m;
        std::vector<Node> nodes;
        std::vector<Link> links;

        std::vector<int> pivots;
        bool deg_rep = false;

        std::vector<std::vector<double>> sp;
        std::vector<constraint> constraints;
        std::vector<constraint> constraintbh;
        //set<constraint_sgd> set_sgd;
        std::vector<std::vector<merge_cons>> append_cons;
        std::vector<std::vector<merge_cons>> append_t_cons;
        bool t_schema=false;
        std::vector<constraint_sgd> constraints_sgd;
        std::vector<force> bh_forces;

        std::vector<t_force> bh_t_forces;
        std::vector<float> count;
        std::vector<float> force_r;
        std::vector<float> schedule;
        std::vector<float> schedulet;
        rk_state rstate;
        para paraSGD;
        std::set<constraint> consBH;
        para paraBH;
        double alphaDecay = 0.034;
        bool using_NP=false;
        double threshold;
    
        int get_n(){return n;}

        void set_n(int _n){n=_n;}

        int get_m(){return m;}

        void set_m(int _m){n=_m;}

        
        graph();

        void readGraph(std::string filename);

        void readGraphMtx(std::string filename);

        void initgraph(std::string filename);
        void initgraph(std::string filename,bool stringid);
        void max_connect();

        void initEdge(int m);

        void initNode(int n);

        void initPosition();

        void initRandomPosition();

        void initPivotMDSPosition(int p = 100);

        void readPosition(std::string filename);

        void readGVPosition(std::string filename);

        void readLinlogPosition(std::string filename);

        void initPosition(std::string filename);

        void solveDij();

        /*   void SMLayout();
           void FDPLayout();
           void FA2Layout();
           void MaxentLayout();
           void BSMLayout();*/

        void append_E_range(para pa);

        void append_N_range(para pa);

        void append_Np_range(para pa);

        void append_S_range(int neighbor, para pa);

        void append_range(std::vector<Link> nodeList, para pa);

        void append_BH_range(int neighbor, para pa);

        void append_E_range(force pa);

        void append_N_range(force pa);

        void append_Np_range(force pa);

        void append_E_range(t_force pa);

        void append_N_range(t_force pa);

        void append_S_range(int neighbor, force pa);

        void append_range(std::vector<Link> nodeList, force pa);

        void append_BH_range2(force pa);

        void append_BH_t_force(t_force pa);

        void preSolve_t_BH(int t_max=200,double decay=0.034f,double threshold=0.0001f);

        void preSolveSGD(double eps = 0.1, int t_max = 200, int seed = 42, float eta_max = 100.0f);

        void preSolveSGD2(double eps = 0.1, int t_max = 200, int seed = 42, float eta_max = 100.0f);
        void preSolveSGD_linear(float e_max=100.0f,float e_min=0.01f, int t_max = 200, int seed = 42);
        void preSolveBH();

        void preSolveBH2();

        bool solveSGD(int iter);

        bool solveSGD2(int iter);

        bool solveSGD_t(int iter);
        bool solve_bilevel(int iter);
        
        bool power_law_graph();

        void solveBH(int iter);

        void solveBH2(int iter);

        void solve_t_BH(int iter);

        void outToFile(std::string filename);

        void RelationOutToFile(std::string filename);

        void outToFile(std::string filename, std::string methodname);

        void RelationOutToFile(std::string filename, std::string methodname);

        void MtxToCsv(std::string filename);

        void drawSVG(std::string filename, std::string method,bool output_pos=false);

    };
#endif //TAURUS_SOLVER_H

   

