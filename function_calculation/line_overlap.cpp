#include "special_functions.h"
#include "integration.h"
#include "constants.h"
#include "line_overlap.h"


line_overlap_func_x2::line_overlap_func_x2() 
    : d(0.), g1(0.), g2(0.), dx(0.), mu(0.), x1(0.), y1(0.), y2(0.), gz1(0.), gz2(0.), dz(0.), lprofile(0)
{;}

void line_overlap_func_x2::set_line_profile(line_profile* lp){
    lprofile = lp;
}

bool line_overlap_func_x2::condition(double t) {
    return (-0.5 * (lprofile->integral(x1 + t) - y1) / gz1 - 0.5 * (lprofile->integral(x1 - dx + t) - y2) / gz2 - t / dz < -50.);
}

// 1/sqrt(pi) is outside the integral:
double line_overlap_func_x2::operator() (double x2) {
    return exp(-x2 * x2 - 0.5 * (lprofile->integral(x2) - y1) / gz1 - 0.5 * (lprofile->integral(x2 - dx) - y2) / gz2 - (x2 - x1) / dz);
}

void line_overlap_func_x2::set_parameters(double dd, double gg1, double gg2, double dxx) { 
    d = dd; 
    g1 = gg1; 
    g2 = gg2; 
    dx = dxx; 
}

void line_overlap_func_x2::set_mu(double m) { 
    mu = m; 
    gz1 = g1 * mu * mu; 
    gz2 = g2 * mu * mu; 
    dz = d * mu * mu; 
}

void line_overlap_func_x2::set_x1(double x) {
    x1 = x;
    y1 = lprofile->integral(x1);
    y2 = lprofile->integral(x1 - dx);
}


line_overlap_func_x1::line_overlap_func_x1(double ta, bool v)
    : tau(ta), verbosity(v), d(0.), b(0.), g1(0.), g2(0.), dx(0.), mu(0.), gz1(0.), gz2(0.), t(0.), 
    lprofile(0), func2(0)
{;}

void line_overlap_func_x1::set_line_profile(line_profile* lp) {
    lprofile = lp;
}

double line_overlap_func_x1::operator() (double x1)
{
    func2->set_x1(x1);
    t = tau * mu; // is auxiliary variable, is not defined by default (tau < 0)

    if (x1 + t > lprofile->x_max || t < DBL_EPSILON)
        t = lprofile->x_max - x1;

    while (func2->condition(t)) { 
        t *= 0.75; 
    }
    return 0.;
}

void line_overlap_func_x1::set_func(line_overlap_func_x2* f) { 
    func2 = f; 
}

void line_overlap_func_x1::set_mu(double m) {
    mu = m;
    gz1 = g1 * mu * mu;
    gz2 = g2 * mu * mu;
    func2->set_mu(mu);
}

void line_overlap_func_x1::set_parameters(double dd, double gg1, double gg2, double dxx) {
    d = dd;
    g1 = gg1;
    g2 = gg2;
    dx = dxx;
    b = g1 / d;
    func2->set_parameters(d, g1, g2, dx);
}


line_overlap_sf1::line_overlap_sf1(double t, bool v) : line_overlap_func_x1(t, v) 
{;}

double line_overlap_sf1::operator() (double x1) 
{
    line_overlap_func_x1::operator()(x1);    
    if (t > 1.e-6) { // 1/sqrt(pi) is outside the integral;
        return (*lprofile)(x1) * qromb<line_overlap_func_x2>(*func2, x1, x1 + t, 1.e-7, verbosity);
    }
    else {
        return lprofile->approx1(x1, dx, gz1, gz2, b);
    }
}


line_overlap_sf2::line_overlap_sf2(double t, bool v) : line_overlap_func_x1(t, v) 
{;}

double line_overlap_sf2::operator()(double x1) 
{
    line_overlap_func_x1::operator()(x1);  
    if (t > 1.e-6) { // 1/sqrt(pi) is outside the integral;
        return (*lprofile)(x1 - dx) * qromb<line_overlap_func_x2>(*func2, x1, x1 + t, 1.e-7, verbosity);
    }
    else {
        return lprofile->approx2(x1, dx, gz1, gz2, b);
    }
}


line_overlap_func_mu::line_overlap_func_mu(bool v) : verbosity(v), d(0.), g1(0.), g2(0.), dx(0.), func1(0) 
{;}

double line_overlap_func_mu::operator()(double mu) 
{
    func1->set_mu(mu);
    return qromb<line_overlap_func_x1>(*func1, func1->lprofile->x_min, func1->lprofile->x_max, 1.e-6, verbosity) / (mu * mu); // gamma is outside the integral;
}

void line_overlap_func_mu::set_func(line_overlap_func_x1* f1) { 
    func1 = f1;  
}

void line_overlap_func_mu::set_parameters(double dd, double gg1, double gg2, double dxx) {
    d = dd;
    g1 = gg1;
    g2 = gg2;
    dx = dxx;
    func1->set_parameters(d, g1, g2, dx);
}
