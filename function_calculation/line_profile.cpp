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
    return exp(-x*x);
}

double line_profile_gauss::integral(double x) {
    return ef.f(x);
}

double line_profile_gauss::approx1(double x1, double dx, double gz1, double gz2, double b)
{
    double e1, e2;
    e1 = (*this)(x1);
    e2 = (*this)(x1 - dx);
    return e1 * e1 * gz1 / (2. * x1 * gz1 + ONEDIVBY_SQRT_PI * (e1 + e2 * gz1 / gz2) + b); // 1/pi is outside the integral;
}

double line_profile_gauss::approx2(double x1, double dx, double gz1, double gz2, double b)
{
    double e1, e2;
    e1 = (*this)(x1);
    e2 = (*this)(x1 - dx);
    return e1 * e2 * gz1 / (2. * x1 * gz1 + ONEDIVBY_SQRT_PI * (e1 + e2 * gz1 / gz2) + b); // 1/pi is outside the integral;
}


line_profile_rectang::line_profile_rectang(double xr) : line_profile(xr)
{
    eq_width = SQRT_PI;
    eq_height = ONEDIVBY_SQRT_PI;
}

double line_profile_rectang::operator() (double x) {
    return (fabs(x) < 0.5 * eq_width) ? SQRT_PI * eq_height : 0.;
}

double line_profile_rectang::integral(double x) {
    if (x < -0.5 * eq_width)
        return 0.;
    else if (x < 0.5 * eq_width)
        return 2. * eq_height * (x + 0.5 * eq_width);
    else 
        return 2.;
}

double line_profile_rectang::approx1(double x1, double dx, double gz1, double gz2, double b)
{
    if (fabs(x1) > 0.5 * eq_width)
        return 0.;
    else {
        double e1, e2;
        e1 = (*this)(x1);
        e2 = (*this)(x1 - dx);
        return e1 * e1 * gz1 / (ONEDIVBY_SQRT_PI * (e1 + e2 * gz1 / gz2) + b);
    }
}

double line_profile_rectang::approx2(double x1, double dx, double gz1, double gz2, double b)
{
    if (fabs(x1 - dx) > 0.5 * eq_width || fabs(x1) > 0.5 * eq_width)
        return 0.;
    else {
        double e1, e2;
        e1 = (*this)(x1);
        e2 = (*this)(x1 - dx);
        return e1 * e2 * gz1 / (ONEDIVBY_SQRT_PI * (e1 + e2 * gz1 / gz2) + b);
    }
}
