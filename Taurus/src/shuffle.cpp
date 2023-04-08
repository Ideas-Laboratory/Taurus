#include "shuffle.h"
long
rk_long(rk_state *state)
{
    return rk_ulong(state) >> 1;
}

unsigned long
rk_ulong(rk_state *state)
{
#if ULONG_MAX <= 0xffffffffUL
    return rk_random(state);
#else
    return (rk_random(state) << 32) | (rk_random(state));
#endif
}
unsigned long
rk_random(rk_state *state)
{

    unsigned long y;

    if (state->pos == RK_STATE_LEN) {
        int i;

        for (i = 0; i < 624 - 397; i++) {
            y = (state->key[i] & 0x80000000UL) | (state->key[i+1] & 0x7fffffffUL);
            state->key[i] = state->key[i+397] ^ (y>>1) ^ (-(y & 1) & 0x9908b0dfUL);
        }
        for (; i < 624 - 1; i++) {
            y = (state->key[i] & 0x80000000UL) | (state->key[i+1] & 0x7fffffffUL);
            state->key[i] = state->key[i+(397-624)] ^ (y>>1) ^ (-(y & 1) & 0x9908b0dfUL);
        }
        y = (state->key[624 - 1] & 0x80000000UL) | (state->key[0] & 0x7fffffffUL);
        state->key[624 - 1] = state->key[397 - 1] ^ (y >> 1) ^ (-(y & 1) & 0x9908b0dfUL);

        state->pos = 0;
    }
    y = state->key[state->pos++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}
unsigned long
rk_interval(unsigned long max, rk_state *state)
{
    unsigned long mask = max, value;

    if (max == 0) {
        return 0;
    }
    /* Smallest bit mask >= max */
    mask |= mask >> 1;
    mask |= mask >> 2;
    mask |= mask >> 4;
    mask |= mask >> 8;
    mask |= mask >> 16;
#if ULONG_MAX > 0xffffffffUL
    mask |= mask >> 32;
#endif

    /* Search a random value in [0..mask] <= max */
#if ULONG_MAX > 0xffffffffUL
    if (max <= 0xffffffffUL) {
        while ((value = (rk_random(state) & mask)) > max);
    }
    else {
        while ((value = (rk_ulong(state) & mask)) > max);
    }
#else
    while ((value = (rk_ulong(state) & mask)) > max);
#endif
    return value;
}
void
rk_seed(unsigned long seed, rk_state *state)
{
    int pos;
    seed &= 0xffffffffUL;

    /* Knuth's PRNG as used in the Mersenne Twister reference implementation */
    for (pos = 0; pos < RK_STATE_LEN; pos++) {
        state->key[pos] = seed;
        seed = (1812433253UL * (seed ^ (seed >> 30)) + pos + 1) & 0xffffffffUL;
    }
    state->pos = RK_STATE_LEN;
    state->gauss = 0;
    state->has_gauss = 0;
    state->has_binomial = 0;
}

void fisheryates_shuffle(std::vector<constraint> &terms, rk_state &rstate)
{
    int n = terms.size();
    //clock_t time_start = clock();
    for (unsigned i=n-1; i>=1; i--)
    {

        unsigned j = rk_interval(i, &rstate);
        constraint temp = terms[i];

        terms[i] = terms[j];
        terms[j] = temp;
        //std::swap(terms[i],terms[j]);
    }

    //clock_t time_end = clock();
    //cout << "rk time:" << 1000 * (time_end - time_start) / (double)CLOCKS_PER_SEC << endl;
}
void fisheryates_shuffle(std::vector<constraint_sgd> &terms, rk_state &rstate)
{
    int n = terms.size();
    //clock_t time_start = clock();
    for (unsigned i=n-1; i>=1; i--)
    {
        unsigned j = rk_interval(i, &rstate);
        constraint_sgd temp = terms[i];

        terms[i] = terms[j];
        terms[j] = temp;
        //std::swap(terms[i],terms[j]);
    }
}
void fisheryates_shuffle(std::vector<int> &terms, rk_state &rstate)
{
    int n = terms.size();
    // clock_t time_start = clock();
    for (unsigned i=n-1; i>=1; i--)
    {
        unsigned j = rk_interval(i, &rstate);
        int temp = terms[i];
        terms[i] = terms[j];
        terms[j] = temp;
        //std::swap(terms[i],terms[j]);
    }
//    clock_t time_end = clock();
//    cout << "rk time:" << 1000 * (time_end - time_start) / (double)CLOCKS_PER_SEC << endl;
}
