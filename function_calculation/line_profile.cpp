#include <cmath>
#include "constants.h"
#include "line_profile.h"

line_profile::line_profile(double xr) : x_min(-xr), x_max(xr)
{;}

void line_profile::set_x_range(double x1, double x2) {
    x_min = x1;
    x_max = x2;
}


line_profile_gauss::line_profile_gauss(double xr) : line_profile(xr)
{;}

double line_profile_gauss::operator() (double x) {
    return ONEDIVBY_SQRT_PI*exp(-x*x);
}

// used in the method, error function (x) = 2/sqrt(pi) *int_{0}^{x} exp(-t*t) dt
double line_profile_gauss::integral(double x) {
    return 0.5*ef.f(x);
}

double line_profile_gauss::approx1(double x1, double dx, double gz1, double gz2, double b)
{
    double e1, e2;
    e1 = (*this)(x1);
    e2 = (*this)(x1 - dx);
    return e1 * e1 * gz1 / (2. * x1 * gz1 + e1 + e2 * gz1 / gz2 + b);
}

double line_profile_gauss::approx2(double x1, double dx, double gz1, double gz2, double b)
{
    double e1, e2;
    e1 = (*this)(x1);
    e2 = (*this)(x1 - dx);
    return e1 * e2 * gz1 / (2. * x1 * gz1 + e1 + e2 * gz1 / gz2 + b);
}

// Rectangular profile
line_profile_rectang::line_profile_rectang() : line_profile(0.)
{
    eq_semiwidth = 0.5*SQRT_PI;
    eq_height = 0.5/eq_semiwidth;
    
    x_min = -eq_semiwidth;
    x_max = eq_semiwidth;
}

double line_profile_rectang::operator() (double x) {
    return (fabs(x) < eq_semiwidth) ? eq_height : 0.;
}

double line_profile_rectang::integral(double x) {
    if (x < -eq_semiwidth)
        return 0.;
    else if (x < eq_semiwidth)
        return eq_height * (x + eq_semiwidth);
    else 
        return 1.;
}

