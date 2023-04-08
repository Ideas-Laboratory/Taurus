#include <random>
#include "RepulsiveForce.h"
FruchtermanReingold::FruchtermanReingold()
{
    grid_quotient(2);
}


void FruchtermanReingold::calculate_exact_repulsive_forces(
        vector<subgraph> &sg,int act_level,vector<DPoint>& F_rep)
{
    //naive algorithm by Fruchterman & Reingold
    DPoint nullpoint(0, 0);
    int node_number = sg[act_level].sub_node.size();
    vector<int> array_of_the_nodes(node_number + 1);

    for (int i = 0; i < node_number; ++i) {
        F_rep[i] = nullpoint;
    }

    int counter = 1;
    for (int i = 0; i < node_number; ++i) {
        array_of_the_nodes[counter] = i;
        counter++;
    }

    for (int i = 1; i < node_number; i++) {
        for (int j = i + 1; j <= node_number; j++)
        {
            int u = array_of_the_nodes[i];
            int v = array_of_the_nodes[j];
            DPoint u_position,v_position;
            u_position.first=sg[act_level].sub_node[u].m_x;
            u_position.second=sg[act_level].sub_node[u].m_y;
            v_position.first=sg[act_level].sub_node[v].m_x;
            v_position.second=sg[act_level].sub_node[v].m_y;
            DPoint f_rep_u_on_v = numexcept::f_rep_u_on_v(u_position, v_position);
            F_rep[v].first += f_rep_u_on_v.first;
            F_rep[v].second += f_rep_u_on_v.second;
            F_rep[u].first -= f_rep_u_on_v.first;
            F_rep[u].second -= f_rep_u_on_v.second;
        }
    }
}


void FruchtermanReingold::make_initialisations(double bl, DPoint d_l_c, int grid_quot)
{
    grid_quotient(grid_quot);
    down_left_corner = d_l_c; //export this two values from FMMM
    boxlength = bl;
}
static std::mt19937 s_random;
int randomNumber(int low, int high)
{

    std::uniform_int_distribution<> dist(low,high);
    return dist(s_random);
}
DPoint numexcept::choose_distinct_random_point_in_disque(DPoint old_point,
                                                         double xmin,double xmax,double ymin,double ymax)
{
    const int BILLION = 1000000000;
    double mindist;//minimal distance from old_point to the boundaries of the disc
    double mindist_to_xmin,mindist_to_xmax,mindist_to_ymin,mindist_to_ymax;
    double rand_x,rand_y;
    DPoint new_point;

    mindist_to_xmin = old_point.first - xmin;
    mindist_to_xmax = xmax -  old_point.first;
    mindist_to_ymin = old_point.second - ymin;
    mindist_to_ymax = ymax -  old_point.second;

    mindist = min(min(mindist_to_xmin,mindist_to_xmax), min(mindist_to_ymin,mindist_to_ymax));

    if (mindist > 0) {
        do {
            //assign random double values in range (-1,1)
            rand_x = 2*(double(randomNumber(1,BILLION)+1)/(BILLION+2)-0.5);
            rand_y = 2*(double(randomNumber(1,BILLION)+1)/(BILLION+2)-0.5);
            new_point.first = old_point.first+mindist*rand_x*epsilon;
            new_point.second = old_point.second+mindist*rand_y*epsilon;
        } while ((old_point == new_point)||(sqrt((old_point.first-new_point.first)*(old_point.first-new_point.first)+(old_point.second-new_point.second)*(old_point.second-new_point.second)) >= mindist*epsilon));
    } else if (mindist == 0) { //old_point lies at the boundaries
        double mindist_x =0;
        double mindist_y =0;

        if (mindist_to_xmin > 0)
            mindist_x = (-1)* mindist_to_xmin;
        else if (mindist_to_xmax > 0)
            mindist_x = mindist_to_xmax;
        if (mindist_to_ymin > 0)
            mindist_y = (-1)* mindist_to_ymin;
        else if (mindist_to_ymax > 0)
            mindist_y = mindist_to_ymax;

        if (mindist_x != 0 || mindist_y != 0) {
            do {
                //assign random double values in range (0,1)
                rand_x = double(randomNumber(1,BILLION)+1)/(BILLION+2);
                rand_y = double(randomNumber(1,BILLION)+1)/(BILLION+2);
                new_point.first = old_point.first+mindist_x*rand_x*epsilon;
                new_point.second = old_point.second+mindist_y*rand_y*epsilon;
            } while (old_point == new_point);
        } else {
            std::cout<<"Error DIM2:: box is equal to old_pos"<<std::endl;
        }
    } else {
        std::cout<<"Error DIM2:: choose_distinct_random_point_in_disque: old_point not ";
        std::cout<<"in box"<<std::endl;
    }

    return new_point;
}

DPoint numexcept::choose_distinct_random_point_in_radius_epsilon(DPoint old_pos)
{
    double xmin = old_pos.first-1*epsilon;
    double xmax = old_pos.first+1*epsilon;
    double ymin = old_pos.second-1*epsilon;
    double ymax = old_pos.second+1*epsilon;

    return choose_distinct_random_point_in_disque(old_pos,xmin,xmax,ymin,ymax);
}

static double random_precision_number(double shift)
{
    const int BILLION = 1000000000;
    double rand = shift + double(randomNumber(1, BILLION) + 1) / (BILLION + 2);
    return randomNumber(0, 1) == 0 ? rand : -rand;
}

bool numexcept::f_rep_near_machine_precision(double distance, DPoint& force)
{
    const double POS_BIG_LIMIT = POS_BIG_DOUBLE * 1e-190;
    const double POS_SMALL_LIMIT = POS_SMALL_DOUBLE * 1e190;

    if (distance > POS_BIG_LIMIT) {
        force = DPoint(POS_SMALL_LIMIT * random_precision_number(1),
                       POS_SMALL_LIMIT * random_precision_number(1));
        return true;
    } else if (distance < POS_SMALL_LIMIT) {
        force = DPoint(POS_BIG_LIMIT * random_precision_number(0),
                       POS_BIG_LIMIT * random_precision_number(0));
        return true;
    }
    return false;
}

bool numexcept::f_near_machine_precision(double distance, DPoint& force)
{
    const double POS_BIG_LIMIT = POS_BIG_DOUBLE * 1e-190;
    const double POS_SMALL_LIMIT = POS_SMALL_DOUBLE * 1e190;

    if (distance < POS_SMALL_LIMIT) {
        force = DPoint(POS_SMALL_LIMIT * random_precision_number(1),
                       POS_SMALL_LIMIT * random_precision_number(1));
        return true;
    } else if (distance > POS_BIG_LIMIT) {
        force = DPoint(POS_BIG_LIMIT * random_precision_number(0),
                       POS_BIG_LIMIT * random_precision_number(0));
        return true;
    }
    return false;
}


bool numexcept::nearly_equal(double a,double b)
{
    double delta = 1e-10;
    double small_b,big_b;

    if(b > 0) {
        small_b = b*(1-delta);
        big_b   = b*(1+delta);

    } else //b <= 0
    {
        small_b = b*(1+delta);
        big_b   = b*(1-delta);
    }

    return (small_b <= a) && (a <= big_b);
}

DPoint numexcept::f_rep_u_on_v(DPoint pos_u, DPoint pos_v)
{
    if (pos_u == pos_v) {
        // Exception handling if two nodes have the same position
        pos_u = choose_distinct_random_point_in_radius_epsilon(pos_u);
    }
    DPoint vector_v_minus_u;
    vector_v_minus_u.first= pos_v.first - pos_u.first;
    vector_v_minus_u.second= pos_v.second - pos_u.second;
    double norm_v_minus_u = sqrt(vector_v_minus_u.first*vector_v_minus_u.first+vector_v_minus_u.second*vector_v_minus_u.second);
    DPoint f_rep_u_on_v;
    if (!f_rep_near_machine_precision(norm_v_minus_u, f_rep_u_on_v)) {
        double scalar = f_rep_scalar(norm_v_minus_u) / norm_v_minus_u;
        f_rep_u_on_v.first = scalar * vector_v_minus_u.first;
        f_rep_u_on_v.second = scalar * vector_v_minus_u.second;
    }
    return f_rep_u_on_v;
}