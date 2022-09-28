

#ifndef TAURUS_LAYOUT_H
#define TAURUS_LAYOUT_H

#endif //TAURUS_LAYOUT_H

#include "Multilevel.h"
void SMLayout(graph &g);
void BSMLayout(graph &g);
void FDPLayout(graph &g);
void FA2Layout(graph &g);
void LinLogLayout(graph &g);
void MaxentLayout(graph &g);

void Layout(graph &g,vector<force> f);
