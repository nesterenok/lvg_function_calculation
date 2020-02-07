#pragma once
#include "line_profile.h"

// J1 = S1 *K11 + S2 * K12
// it is assumed that dx = (frequency_1 - frequency_2) / d_freq_doppler, difference between frequency_1 and frequency_2 (except dx) is omitted,
// gz1 and gz2 (or velocity gradient) are positive, lower opacity - larger g or d, 
class line_overlap_func2
{
public:
    line_profile* lprofile;
    double d, g1, g2, dx, mu, x1, y1, y2, gz1, gz2, dz;

    // the tabulated (external) parameters: 
    void set_parameters(double dd, double gg1, double gg2, double dxx);
    
    // the integration variables:
    void set_mu(double m);
    void set_x1(double x);
    void set_line_profile(line_profile* lp);

    virtual double operator() (double x2);

    // true - the value of t is small enough
    virtual bool condition(double t);
    line_overlap_func2();
};

// basic class for integrated functions:
class line_overlap_func1
{
public:
    bool verbosity;
    double d, b, g1, g2, dx, mu, gz1, gz2, t, tau;
    
    line_profile* lprofile;
    line_overlap_func2* func2;
   
    void set_line_profile(line_profile* lp);
    void set_func(line_overlap_func2* f);
    void set_parameters(double dd, double gg1, double gg2, double dxx);
    void set_mu(double m);

    virtual double operator() (double x1);

    // tau is gamma*k_L*z = dv/dz *z/v_th, z is length of the region in question;
    line_overlap_func1(double ta, bool v);
};

// the first integral on x2
// 1./sqrt(pi) exp(-x1*x1) *\int_{x1}^{x1+tau*mu} dx2 func2(x2)
class line_overlap_sf1 : public line_overlap_func1
{
public:
    double operator() (double x1);
    line_overlap_sf1(double t = -1, bool v = false);
};

class line_overlap_sf2 : public line_overlap_func1
{
public:
    double operator() (double x1);
    line_overlap_sf2(double t = -1, bool v = false);
};

// Integral on x1, 1/(mu*mu) *\int_{-\infty}^{\infty} dx1 func1(x1)
class line_overlap_func_mu
{
public:
    bool verbosity;
    double d, g1, g2, dx, x1_min, x1_max; // integration limits for external integration on x, the variable can be shifted;
    
    line_overlap_func1* func1;

    void set_func(line_overlap_func1* f1);
    void set_parameters(double dd, double gg1, double gg2, double dxx);
    void set_integr_lim(double x1_min, double x1_max);
    double operator() (double mu);

    line_overlap_func_mu(bool v = false);
};

// There is distinct class for rectangular profile,
class line_overlap_rect 
{
public:
    double d, g1, g2, dx;
    line_profile_rectang* lprofile;
    
    void set_parameters(double dd, double gg1, double gg2, double dxx);
    line_overlap_rect();
};

class line_overlap_rect1 : public line_overlap_rect
{
public: 
    double operator() (double mu);
};

class line_overlap_rect2 : public line_overlap_rect
{
public:
    double operator() (double mu);
};

