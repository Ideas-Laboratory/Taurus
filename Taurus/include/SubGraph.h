//
// Created by xml on 2022/8/11.
//

#ifndef TEST_CPP_SUBGRAPH_H
#define TEST_CPP_SUBGRAPH_H

#include "vector"

class subNode{
public:
    float m_x,m_y;
    double weight,height;
    subNode* v_lower_level;
    subNode* v_higher_level;
    int v_higher_level_index,v_lower_level_index;
    int id;
    int mass;
    int type;
    int dedicated_sun_node;//同层的sun_node
    double dedicated_sun_distance;
    int dedicated_pm_node;
    std::vector<double> lambda; //!< the factors lambda for scaling the length of this edge
    //! relative to the pass between v's sun and the sun of a
    //! neighbour solar system

    std::vector<int> neighbour_s_node;//!< this is the list of the neighbour solar systems suns
    //! lambda[i] corresponds to neighbour_s_node[i]

    //vector<double>* lambda_List_ptr; //!< a poi
    //vector<int>* neighbour_s_node_List_ptr; //!< a pointer to to the neighbour_s_node list
    std::vector<int>  moon_List;//!< the list of all dedicated moon nodes (!= nil if type == 3)
    //vector<int>* moon_List_ptr;//!< a pointer to the moon_List
    bool placed;   //!< indicates weather an initial position has been assigned to this
    //! node or not
    double angle_1;//!< describes the sector where nodes that are not adjacent to other
    double angle_2;//!< solar systems have to be placed

};
class adjEdge{
public:
    int st,end;
    int st_index,end_index;
    bool moon_edge=false;
    double length;
};
class adjLink{
public:
    int source,target;
    mutable double length;
    mutable int count=0;
    adjLink* e_original;
    adjLink* e_subgraph;
    adjLink() {}

    adjLink(int s, int t) : source(s), target(t) {}
    adjLink(int s, int t,double l,int c) : source(s), target(t) ,length(l),count(c){}

    bool operator<(const adjLink &al) const {
        if (source!= al.source) {
            return source < al.source;
        } else {
            return target < al.target;
        }
    }
};
class subgraph{
public:
    int numofnode,numofedge;
    std::vector<subNode> sub_node;
    std::vector<adjLink> sub_edge;
    std::vector<std::vector<adjEdge>> edge;
};
#endif //TEST_CPP_SUBGRAPH_H