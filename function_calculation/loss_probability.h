#pragma once
#include "special_functions.h"
#include "integration.h"
#include "constants.h"

#define X_RANGE 7. 

// used in the method, error function (x) = 2/sqrt(pi) *int_{0}^{x} exp(-t*t) dt
class lvg_func_x2
{
public:
    double d, g, mu, x1, y1, gz, dz;
    err_func ef;

    void set_parameters(double dd, double gg) { d = dd; g = gg; };
    void set_mu(double m) { mu = m; gz = g * mu * mu; dz = d * mu * mu; };
    void set_x1(double x) { x1 = x; y1 = ef.f(x1); };

    // 1/sqrt(pi) is outside the integral:
    virtual double operator() (double x2) {
        return exp(-x2 * x2 - 0.5 * (ef.f(x2) - y1) / gz - (x2 - x1) / dz);
    };
    // true - the value of t is small enough
    virtual bool condition(double t) {
        return (-0.5 * (ef.f(x1 + t) - y1) / gz - t / dz < -50.);
    };
};

// basic class for integrated functions:
class lvg_func_x1
{
public:
    double b, d, g, e, mu, gz, t, tau, upper_lim;
    lvg_func_x2* func2;
    err_func ef;

    void set_func(lvg_func_x2* f) { func2 = f; };
    void set_parameters(double dd, double gg) {
        d = dd;
        g = gg;
        b = g / d;
        func2->set_parameters(d, g);
    };
    void set_mu(double m) {
        mu = m;
        gz = g * mu * mu;
        func2->set_mu(mu);
    };
    virtual double operator() (double x1)
    {
        func2->set_x1(x1);
        t = tau * mu;

        if (x1 + t > upper_lim || t < 0.)
            t = upper_lim - x1;

        while (func2->condition(t)) { t *= 0.75; }
        return 0.;
    }
    // tau is gamma*k_L*z = dv/dz *z/v_th, z is length of the region in question,
    lvg_func_x1(double t) : b(0.), d(0.), g(0.), e(0.), mu(0.), gz(0.), t(0.), tau(t), upper_lim(X_RANGE) { ; } 
};

// function to calculate one-sided loss probability for photons created by line processes:
// exp(-x1*x1) *\int_{x1}^{x1+tau*mu} dx2 func2(x2)
class lvg_func_line_x1 : public lvg_func_x1
{
public:
    double operator() (double x1) {
        lvg_func_x1::operator()(x1);

        if (t > 1.e-6) // 1/sqrt(pi) is outside the integral;
            return exp(-x1 * x1) * qromb<lvg_func_x2>(*func2, x1, x1 + t, 1.e-7);
        else {
            e = exp(-x1 * x1);
            return e * e * gz / (2. * x1 * gz + ONEDIVBY_SQRT_PI * e + b); // 1/pi is outside the integral;
        }
    };
    lvg_func_line_x1(double t = -1) : lvg_func_x1(t) { ; }
};

// function to calculate one-sided loss probability for photons created by continuum processes:
// \int_{x1}^{x1+tau*mu} dx2 func2(x2)
class lvg_func_cont_x1 : public lvg_func_x1
{
public:
    double operator() (double x1)
    {
        lvg_func_x1::operator()(x1);

        if (t > 1.e-6)
            return qromb<lvg_func_x2>(*func2, x1, x1 + t, 1.e-7);
        else {
            e = exp(-x1 * x1);
            return gz * e / (2. * x1 * gz + ONEDIVBY_SQRT_PI * e + b); // 1/sqrt(pi) is outside the integral;
        }
    };
    lvg_func_cont_x1(double t = -1) : lvg_func_x1(t) { ; }
};

// Integral on x1, 1/(mu*mu) *\int_{-\infty}^{\infty} dx1 func1(x1)
struct lvg_func_mu
{
    double d, g;
    lvg_func_x1* func1;

    void set_func(lvg_func_x1* f1, lvg_func_x2* f2) { func1 = f1; func1->set_func(f2); };

    void set_parameters(double dd, double gg) {
        d = dd;
        g = gg;
        func1->set_parameters(d, g);
    };
    double operator() (double mu) {
        func1->set_mu(mu);
        return qromb<lvg_func_x1>(*func1, -X_RANGE, X_RANGE, 1.e-6) / (mu * mu); // gamma is outside the integral;
    };
};

// expression (2.20) Hummer, Rybicki ApJ 293, p.258 (1985);
// loss probability function for photon in the absence of continuum absorption (must be multiplied by 0.5 to get one-sided);
class lvg_func_p0 {
public:
    double g;
    double operator() (double mu) {
        return mu * mu * (1. - exp(-1. / (g * mu * mu)));
    };
};

//
class lvg_aux {
public:
    double b;
    double operator() (double x) {
        return ONEDIVBY_SQRT_PI * exp(-x * x) / (ONEDIVBY_SQRT_PI * exp(-x * x) + b);
    };
};

