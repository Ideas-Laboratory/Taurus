#ifndef _RANGE_H_
#define _RANGE_H_
#include "Force.h"
#include "vector"
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

    struct constraint_sgd {
        int i, j;
        double d, w;
        mutable std::vector<force> force_sgd;
        mutable std::vector<t_force> force_t;
        constraint_sgd() {}

        constraint_sgd(int _i, int _j) : i(_i), j(_j), d(0), w(0), force_sgd() {}

        constraint_sgd(int _i, int _j, double _d, double _w, std::vector<force> _force_sgd) : i(_i), j(_j), d(_d), w(_w),
                                                                                         force_sgd(_force_sgd) {}

        constraint_sgd(int _i, int _j, double _d, double _w, std::vector<t_force> _force_t) : i(_i), j(_j), d(_d), w(_w),force_t(_force_t) {}
        bool operator<(const constraint_sgd &_cs) const {
            if (i != _cs.i) {
                return i < _cs.i;
            } else {
                return j < _cs.j;
            }
        }
    };

    struct merge_cons {
        std::vector<force> force_m;
        std::vector<t_force> t_force_m;
        bool flag = false;
        double d, w;
        merge_cons():flag(false){};
        merge_cons(double _d,double _w,std::vector<force> _force):d(_d),w(_w),force_m(_force){}
        merge_cons(double _d,double _w,std::vector<t_force> _force):d(_d),w(_w),t_force_m(_force){}
    };
    struct merge_t_cons {
        std::vector<t_force> force_m;
        bool flag = false;
        double d, w;
        merge_t_cons():flag(false){};
        merge_t_cons(double _d,double _w,std::vector<t_force> _force):d(_d),w(_w),force_m(_force){}
    };
    #endif