#pragma once
#include "special_functions.h"

class line_profile
{
public:
    double x_min, x_max;
    virtual double operator() (double x) = 0; // = normalized profile(x)
    virtual double integral(double x) = 0;    // = integrated normalized profile + const(arbitrary)

    // gz1 = mu*mu*g1, gz2 = mu*mu*g2, b = g1/d
    virtual double approx1(double x1, double dx, double gz1, double gz2, double b) { return 0.; };
    virtual double approx2(double x1, double dx, double gz1, double gz2, double b) { return 0.; };
    
    void set_x_range(double x_min, double x_max);
    line_profile(double x_range);
};

class line_profile_gauss : public line_profile
{
public:
    err_func ef;
    double operator() (double x);
    double integral(double x);
    
    double approx1(double x1, double dx, double gz1, double gz2, double b);
    double approx2(double x1, double dx, double gz1, double gz2, double b);
    
    line_profile_gauss(double x_range = 0.);
};

// Rectangular profile
class line_profile_rectang : public line_profile
{
public:
    double eq_semiwidth, eq_height;
    double operator() (double x);
    double integral(double x);
    line_profile_rectang();
};
