#ifndef _SHUFFLE_H_
#define _SHUFFLE_H_
#include "randomkit.h"
#include "vector"
#include "range.h"
long
rk_long(rk_state *state);

unsigned long
rk_ulong(rk_state *state);

unsigned long
rk_random(rk_state *state);

unsigned long
rk_interval(unsigned long max, rk_state *state);

void
rk_seed(unsigned long seed, rk_state *state);
void fisheryates_shuffle(std::vector<constraint> &terms, rk_state &rstate);
void fisheryates_shuffle(std::vector<constraint_sgd> &terms, rk_state &rstate);
void fisheryates_shuffle(std::vector<int> &terms, rk_state &rstate);
#endif