#include "special_functions.h"
#include "integration.h"
#include "constants.h"
#include "line_overlap.h"


line_overlap_func2::line_overlap_func2() 
    : d(0.), g1(0.), g2(0.), dx(0.), mu(0.), x1(0.), y1(0.), y2(0.), gz1(0.), gz2(0.), dz(0.), lprofile(0)
{;}

void line_overlap_func2::set_line_profile(line_profile* lp){
    lprofile = lp;
}

bool line_overlap_func2::condition(double t) {
    return (-(lprofile->integral(x1 + t) - y1) / gz1 - (lprofile->integral(x1 - dx + t) - y2) / gz2 - t / dz < -50.);
}

// 1/sqrt(pi) is in the integral:
double line_overlap_func2::operator() (double x2) {
    return (*lprofile)(x2) * exp(-(lprofile->integral(x2) - y1) / gz1 - (lprofile->integral(x2 - dx) - y2) / gz2 - (x2 - x1) / dz);
}

void line_overlap_func2::set_parameters(double dd, double gg1, double gg2, double dxx) { 
    d = dd; 
    g1 = gg1; 
    g2 = gg2; 
    dx = dxx; 
}

void line_overlap_func2::set_mu(double m) { 
    mu = m; 
    gz1 = g1 * mu * mu; 
    gz2 = g2 * mu * mu; 
    dz = d * mu * mu; 
}

void line_overlap_func2::set_x1(double x) {
    x1 = x;
    y1 = lprofile->integral(x1);
    y2 = lprofile->integral(x1 - dx);
}


line_overlap_func1::line_overlap_func1(double ta, bool v)
    : tau(ta), verbosity(v), d(0.), b(0.), g1(0.), g2(0.), dx(0.), mu(0.), gz1(0.), gz2(0.), t(0.), 
    lprofile(0), func2(0)
{;}

void line_overlap_func1::set_line_profile(line_profile* lp) {
    lprofile = lp;
}

double line_overlap_func1::operator() (double x1)
{
    func2->set_x1(x1);
    t = tau * mu; // is auxiliary variable, is not defined by default (tau < 0)

    if (x1 + t > lprofile->x_max || t < DBL_EPSILON)
        t = lprofile->x_max - x1;
    
    if (t < 0.)
        t = 0.;

    while (func2->condition(t)) { 
        t *= 0.75; 
    }
    return 0.;
}

void line_overlap_func1::set_func(line_overlap_func2* f) { 
    func2 = f; 
}

void line_overlap_func1::set_mu(double m) {
    mu = m;
    gz1 = g1 * mu * mu;
    gz2 = g2 * mu * mu;
    func2->set_mu(mu);
}

void line_overlap_func1::set_parameters(double dd, double gg1, double gg2, double dxx) {
    d = dd;
    g1 = gg1;
    g2 = gg2;
    dx = dxx;
    b = g1 / d;
    func2->set_parameters(d, g1, g2, dx);
}


line_overlap_sf1::line_overlap_sf1(double t, bool v) : line_overlap_func1(t, v) 
{;}

double line_overlap_sf1::operator() (double x1) 
{
    line_overlap_func1::operator()(x1);    
    if (t > 1.e-6) { // 1/sqrt(pi) is in the integral;
        return (*lprofile)(x1) * qromb<line_overlap_func2>(*func2, x1, x1 + t, 1.e-7, verbosity);
    }
    else {
        return lprofile->approx1(x1, dx, gz1, gz2, b);
    }
}


line_overlap_sf2::line_overlap_sf2(double t, bool v) : line_overlap_func1(t, v) 
{;}

double line_overlap_sf2::operator()(double x1) 
{
    line_overlap_func1::operator()(x1);  
    if (t > 1.e-6) { // 1/sqrt(pi) is in the integral;
        return (*lprofile)(x1 - dx) * qromb<line_overlap_func2>(*func2, x1, x1 + t, 1.e-7, verbosity);
    }
    else {
        return lprofile->approx2(x1, dx, gz1, gz2, b);
    }
}


line_overlap_func_mu::line_overlap_func_mu(bool v) : verbosity(v), d(0.), g1(0.), g2(0.), dx(0.), func1(0), x1_min(0.), x1_max(0.)
{;}

double line_overlap_func_mu::operator()(double mu) 
{
    func1->set_mu(mu);
    return qromb<line_overlap_func1>(*func1, x1_min, x1_max, 1.e-6, verbosity) / (mu * mu); // gamma is outside the integral;
}

void line_overlap_func_mu::set_func(line_overlap_func1* f1) { 
    func1 = f1;  
}

void line_overlap_func_mu::set_parameters(double dd, double gg1, double gg2, double dxx) {
    d = dd;
    g1 = gg1;
    g2 = gg2;
    dx = dxx;
    func1->set_parameters(d, g1, g2, dx);
}

void line_overlap_func_mu::set_integr_lim(double a, double b) {
    x1_min = a;
    x1_max = b;
}

// Rectangular profile
line_overlap_rect::line_overlap_rect() : d(0.), g1(0.), g2(0.), dx(0.) {
    lprofile = new line_profile_rectang();
}


void line_overlap_rect::set_parameters(double dd, double gg1, double gg2, double dxx)
{
    d = dd;  
    g1 = gg1;
    g2 = gg2;
    dx = dxx; 
}

double line_overlap_rect1::operator()(double mu)
{
    double a, b, u, v, dxx, z1, z2, zd, answ(0.);
    
    dxx = fabs(dx); // the answer does not depend on the sign
    z1 = 1./(g1 * mu * mu);
    zd = 1./(d * mu * mu);
    a = lprofile->eq_height * z1 + zd;

    if (dxx < 2. * lprofile->eq_semiwidth) {  
        z2 = 1./(g2 * mu * mu);
        b = a + lprofile->eq_height * z2;
        
        u = 1. - exp(-a * dxx);
        v = 1. - exp(-b * (2. * lprofile->eq_semiwidth - dxx));

        answ = -u / (a * a) + dxx / a;
        answ += u * v / (a * b);
        answ += -v / (b * b) + (2. * lprofile->eq_semiwidth - dxx) / b;
    }
    else {
        answ = -(1. - exp(-a * 2. * lprofile->eq_semiwidth)) / (a * a) + 2. * lprofile->eq_semiwidth / a;
    }  
    return answ * lprofile->eq_height * lprofile->eq_height / (mu * mu);
}

double line_overlap_rect2::operator()(double mu)
{
    double a, b, c, u, v, w, z1, z2, zd, answ(0.);

    z1 = 1. / (g1 * mu * mu);
    z2 = 1. / (g2 * mu * mu);
    zd = 1. / (d * mu * mu);

    if (dx >= 0.) {
        if (dx < 2. * lprofile->eq_semiwidth)
        {
            b = lprofile->eq_height * (z1 + z2) + zd;
            v = 1. - exp(-b * (2. * lprofile->eq_semiwidth - dx));
            answ = -v / (b * b) + (2. * lprofile->eq_semiwidth - dx) / b;
        }
        else answ = 0.;
    }
    else {
        if (dx > -2. * lprofile->eq_semiwidth) {
            a = lprofile->eq_height * z1 + zd;
            c = lprofile->eq_height * z2 + zd;
            b = a + lprofile->eq_height * z2;

            w = 1. - exp(a * dx);
            u = 1. - exp(c * dx);
            v = 1. - exp(-b * (2. * lprofile->eq_semiwidth + dx));
            
            answ = u * v / (b * c);
            answ += -v/(b*b) + (2. * lprofile->eq_semiwidth + dx)/b;
            answ += u * w * (1. - v) / (a * c);
            answ += w * v / (a * b);
        }
        else {
            a = lprofile->eq_height * z1 + zd;
            c = lprofile->eq_height * z2 + zd;
            answ = (1. - exp(-a* 2. * lprofile->eq_semiwidth)) * (1. - exp(-c* 2. * lprofile->eq_semiwidth))
                *exp(zd * (2. * lprofile->eq_semiwidth + dx)) / (a * c);
        }
    }
    return answ * lprofile->eq_height * lprofile->eq_height / (mu * mu);
}
