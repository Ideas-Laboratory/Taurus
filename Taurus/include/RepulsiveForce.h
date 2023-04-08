
#ifndef TEST_CPP_REPULSIVEFORCE_H
#define TEST_CPP_REPULSIVEFORCE_H


#include "SubGraph.h"
#include "vector"
#include "math.h"
#include "iostream"
#define epsilon 0.1
#define POS_SMALL_DOUBLE 1e-300
#define POS_BIG_DOUBLE   1e+300
using namespace std;
typedef pair<float,float> DPoint;
class  FruchtermanReingold
        {
                public:
                //! Constructor
                FruchtermanReingold();

                //! Calculate exact rep. forces for each node.
                void calculate_exact_repulsive_forces(vector<subgraph> &sg,int act_level,vector<DPoint>& F_rep);

                //! Grid approximation of rep.forces for each node.
                void calculate_approx_repulsive_forces(
                        vector<subgraph> &sg,int act_level,vector<DPoint>& F_rep);

                //! Make all initialisations that are needed for FruchtermanReingold.
                void make_initialisations (
                double boxlength,
                DPoint down_left_corner,
                int grid_quotient);

                //! Import updated information of the drawing area.
                void update_boxlength_and_cornercoordinate(double b_l, DPoint d_l_c) {
                    boxlength = b_l; down_left_corner = d_l_c;
                }

                private:
                int _grid_quotient;//!< for coarsening the FrRe-grid
                int max_gridindex; //!< maximum index of a grid row/column
                double boxlength;  //!< length of drawing box
                DPoint down_left_corner;//!< down left corner of drawing box

                //! The number k of rows and colums of the grid is sqrt(|V|) / frGridQuotient()
                //! (Note that in [FrRe] frGridQuotient() is 2.)
                void grid_quotient(int p) { _grid_quotient = ((0<=p) ? p : 2);}
                int grid_quotient() const {return _grid_quotient;}
        };
class numexcept
{
public:
    //! Returns a distinct random point within the smallest disque D with center
    //! old_point that is contained in the box defined by xmin,...,ymax; The size of
    //! D is shrunk by multiplying with epsilon = 0.1; Precondition:
    //! old_point is contained in the box and the box is not equal to old_point.
    static DPoint choose_distinct_random_point_in_disque(
            DPoint old_point,
            double xmin,
            double xmax,
            double ymin,
            double ymax);

    //! If distance has a value near the machine precision the (attractive)force
    //! calculation is not possible (calculated values exceed the machine accuracy) in
    //! this cases true is returned and force is set to a reasonable value that does
    //! not cause problems; Else false is returned and force keeps unchanged.
    static bool f_near_machine_precision(double distance, DPoint& force);

    //! Returns true if a is "nearly" equal to b (needed, when machine accuracy is
    //! insufficient in functions well_seperated and bordering of NMM)
    static bool nearly_equal(double a, double b);

    static DPoint f_rep_u_on_v(DPoint pos_u, DPoint pos_v);

protected:
    //! A random point (distinct from old_pos) on the disque around old_pos with
    //! radius epsilon = 0.1 is computed.
    static DPoint choose_distinct_random_point_in_radius_epsilon(DPoint old_pos);

    //! If distance has a value near the machine precision the repulsive force calculation
    //! is not possible (calculated values exceed the machine accuracy) in this cases
    //! true is returned and force is set to a reasonable value that does
    //! not cause problems; Else false is returned and force keeps unchanged.
    static bool f_rep_near_machine_precision(double distance, DPoint& force);

    //! Returns the repulsing force_function_value of scalar d.
    static double f_rep_scalar(double d) {
        //static_assert(d != 0);
        if(d==0) return 0;
        return 1/d;
    }
};
#endif //TEST_CPP_REPULSIVEFORCE_H