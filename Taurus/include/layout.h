

#ifndef TAURUS_LAYOUT_H
#define TAURUS_LAYOUT_H



#include "Multilevel.h"
void SMLayout(graph &g);
void BSMLayout(graph &g);
void FDPLayout(graph &g);
void FA2Layout(graph &g);
void LinLogLayout(graph &g);
void MaxentLayout(graph &g);
void MaxentLayout(graph &g,int k_neigh);
void Layout(graph &g,vector<force> f);
void t_Layout(graph &g,vector<t_force> f);
#endif 
