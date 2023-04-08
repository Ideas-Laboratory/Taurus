 #ifndef _FORCE_H_
 #define _FORCE_H_
 struct force {
        float range, k, a, b;

        force() {}

        force(float _k, float _a, float _b) : range(-1), k(_k), a(_a), b(_b) {}

        force(float _range, float _k, float _a, float _b) : range(_range), k(_k), a(_a), b(_b) {}
    };

    struct t_force {
        float range, k, a, b, t;

        t_force() {}

        t_force(float _k, float _a, float _b,float _t) : range(-1), k(_k), a(_a), b(_b), t(_t){}

        t_force(float _range, float _k, float _a, float _b,float _t) : range(_range), k(_k), a(_a), b(_b) ,t(_t){}
    };
 struct para {
        float ka, q, a, kr, r, p, pmds;
        double w1 = 1, w2 = 1, w3 = 1;

        para() : ka(0), q(0), a(0), kr(0), r(0), p(0), pmds(0) {}

        para(float ka, float q, float a, float kr, float r, float p, float neighbor, double w1, double w2, double w3)
                : ka(ka), q(q), a(a), kr(kr), r(r), p(p), pmds(neighbor), w1(w1), w2(w2), w3(w3) {}

        para(float ka, float q, float a, float kr, float r, float p, float neighbor) : ka(ka), q(q), a(a), kr(kr), r(r),
                                                                                       p(p), pmds(neighbor) {}
    };
    #endif