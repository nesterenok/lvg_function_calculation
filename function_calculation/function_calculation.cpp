// 
// 31.03.2017. The error was found in approximation formula.
// 05.09.2017. Check for errors. 
//		It takes about 2*24 hours to calculate all p and q integrals using 32 processors, relative accuracy 1.e-5, integration on mu is [1.e-6, 1];
//		Calculation time strongly depends on integration limits for mu;
// 17.03.2020. Check for errors.

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

// the values are used between 1.e-8 and 1.e-2
// the photons with the mu < mu_min must be taken into account
#define MU_MIN 1.e-6        
#define X_RANGE_LINE_OVERLAP 10.

#ifdef __linux__
#include <sys/types.h>
#include <signal.h>
#include <unistd.h>
#elif defined _WIN32 || defined _WIN64 
#include "stdafx.h"
#include <process.h>
#endif

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <cfloat>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <ctime>

#include "utils.h"
#include "integration.h"
#include "loss_probability.h"
#include "line_overlap.h"
#include "line_profile.h"

using namespace std;

void calc_loss_func_p(double g_min, double g_max, double d_min, double d_max, int point_nb_per_order);
void calc_loss_func_q(double g_min, double g_max, double d_min, double d_max, int point_nb_per_order);

// line overlap,
// g_min must be less than 1.e-2,
void calc_line_overlap1(double d_min, double d_max, double g_min, double g_max, double gratio_min, double gratio_max, double dx_min, double dx_max, 
    int point_nb_per_order, int nb_dx_points, bool rect_profile);
void calc_line_overlap2(double d_min, double d_max, double g_min, double g_max, double gratio_min, double gratio_max, double dx_min, double dx_max, 
    int point_nb_per_order, int nb_dx_points, bool rect_profile);

// both min and max values must be positive,
int init_log_grid(double*& arr, double min_val, double max_val, double nb_points_bin)
{
    int i, nb_ep;
    double factor;

    factor = pow(10., 1. / nb_points_bin);
    nb_ep = (int)(log10(max_val / min_val) * nb_points_bin) + 1;

    arr = new double[nb_ep];
    arr[0] = min_val;

    for (i = 1; i < nb_ep; i++) {
        arr[i] = arr[i - 1] * factor;
    }
    return nb_ep;
}

void init_lin_grid(double*& arr, double min_val, double max_val, int nb_points)
{
    int i;
    double factor;

    factor = (max_val - min_val)/ ((double) nb_points - 1.);
    arr = new double[nb_points];
    arr[0] = min_val;

    for (i = 1; i < nb_points; i++) {
        arr[i] = arr[i - 1] + factor;
    }
}

int main(int argc, char** argv)
{

#ifdef _OPENMP
#ifdef __linux__
	omp_set_num_threads(40);
#else
	omp_set_num_threads(1);
#endif

#pragma omp parallel 
	{
#pragma omp master 
		{
			cout << "OpenMP is supported" << endl;
			cout << "Nb of threads: " << omp_get_num_threads() << endl;
		}
	}
#endif

    bool rect_profile;
    int point_nb_per_order, nb_dx_points;
    double g_min, g_max, d_min, d_max, gratio_min, gratio_max, dx_min, dx_max;
//	calc_loss_func_p(g_min = 1.e-6, g_max = 1.e+14, d_min = 1.e-2, d_max = 1.e+9, point_nb_per_order = 16);
//	calc_loss_func_q(g_min = 1.e-6, g_max = 1.e+14, d_min = 1.e-2, d_max = 1.e+9, point_nb_per_order = 16);
    calc_line_overlap1(d_min = 1.e+2, d_max = 1.e+6, g_min = 1.e-6, g_max = 1.e+6, gratio_min = 0.01, gratio_max = 100., dx_min = -4., dx_max = 4., 
        point_nb_per_order = 8, nb_dx_points = 33, rect_profile = false);
    
    calc_line_overlap2(d_min = 1.e+2, d_max = 1.e+6, g_min = 1.e-6, g_max = 1.e+6, gratio_min = 0.01, gratio_max = 100., dx_min = -4., dx_max = 4., 
        point_nb_per_order = 8, nb_dx_points = 33, rect_profile = false);

	return 0;
}


void calc_loss_func_p(double g_min, double g_max, double d_min, double d_max, int point_nb_per_order)
{
	int i, j, nb_g, nb_d;
	double *gamma_arr, *delta_arr;
	
    time_t timer;
	stringstream ss;
	ofstream out1, out2;
	
#ifdef __linux__
	ss.str("");
	ss << "./output/out_p.txt";

	std::ofstream outerr;
	std::streambuf *orig_cerr;
	
	std::ofstream out;
	std::streambuf *orig_cout; 

	outerr.open(ss.str().c_str(), std::ios::app);		// "/dev/null"
	orig_cerr = std::cerr.rdbuf(outerr.rdbuf());	// "/dev/null"
	
	out.open(ss.str().c_str(), std::ios::app);
	orig_cout = std::cout.rdbuf(out.rdbuf());
#endif

	cout << scientific;
	cout.precision(6);

	nb_g = init_log_grid(gamma_arr, g_min, g_max, point_nb_per_order);
	nb_d = init_log_grid(delta_arr, d_min, d_max, point_nb_per_order);

    ss.str("");
    ss << "./output/lvg_loss_func_p.txt";

    out1.open(ss.str().c_str(), std::ios_base::out);
    out1 << scientific;
    out1.precision(5);

    out1 << left << setw(18) << "";
    for (j = 0; j < nb_g; j++) {
        out1 << left << setw(13) << gamma_arr[j];
    }
    out1 << endl;

    ss.str("");
    ss << "./output/lvg_loss_func_p(b=0).txt";

    out2.open(ss.str().c_str(), std::ios_base::out);
    out2 << scientific;

	// Calculation of the loss probability functions;
    timer = time(NULL);
#pragma omp parallel private(i, j)
	{
		lvg_func_p0 calc_p0;
		lvg_func_mu calc_p;

		calc_p.set_func(new lvg_func_line_x1(), new lvg_func_x2());
		
		double *p0 = new double [nb_g];
		memset(p0, 0, nb_g*sizeof(double));

		double **p1 = alloc_2d_array<double>(nb_d, nb_g);
		memset(*p1, 0, nb_d*nb_g*sizeof(double));

#pragma omp for schedule(dynamic, 1), ordered
		for (i = 0; i < nb_d; i++) {
            for (j = 0; j < nb_g; j++)
		    {			
				calc_p.set_parameters(delta_arr[i], gamma_arr[j]);
				// one sided loss probability function for photon created in line processes is 0.5-0.5*p, where p is:
				p1[i][j] = qromb<lvg_func_mu>(calc_p, MU_MIN, 1., 1.e-5, false) /(M_PI*gamma_arr[j]);
			}
			
			#pragma omp ordered
			{
				cout <<left << setw(5) << omp_get_thread_num() << setw(15) << delta_arr[i] 
                    << "time (s): " << (int) (time(NULL)-timer) << endl;

                out1 << left << setw(5) << i << setw(13) << delta_arr[i];
                for (j = 0; j < nb_g; j++) {
                    out1 << left << setw(13) << p1[i][j];
                }
                out1 << endl;
			}
		}
		
#pragma omp for schedule(dynamic, 1), ordered
        for (j = 0; j < nb_g; j++)
        {
            calc_p0.g = gamma_arr[j];
            p0[j] = gamma_arr[j] * qromb<lvg_func_p0>(calc_p0, MU_MIN, 1., 1.e-8);

#pragma omp ordered
            {
                out2.precision(6);
                out2 << left << setw(15) << gamma_arr[j] << setw(18) << setprecision(9) << p0[j] << endl; // be carefull: one has to compare 1.-p and p0;
            }
        }
		delete [] p0;
		free_2d_array(p1);
	}
	out1.close();
    out2.close();

	delete [] gamma_arr;
	delete [] delta_arr;
}

void calc_loss_func_q(double g_min, double g_max, double d_min, double d_max, int point_nb_per_order)
{
	int i, j, nb_g, nb_d;
	double *gamma_arr, *delta_arr;

    stringstream ss;
	ofstream output;
		
#ifdef __linux__
	ss.str("");
	ss << "./output/out_q.txt";

	std::ofstream outerr;
	std::streambuf *orig_cerr;
	
	std::ofstream out;
	std::streambuf *orig_cout; 

	outerr.open(ss.str().c_str(), std::ios::app);		// "/dev/null"
	orig_cerr = std::cerr.rdbuf(outerr.rdbuf());	// "/dev/null"
	
	out.open(ss.str().c_str(), std::ios::app);
	orig_cout = std::cout.rdbuf(out.rdbuf());
#endif

	cout << scientific;
	cout.precision(6);

	nb_g = init_log_grid(gamma_arr, g_min, g_max, point_nb_per_order);		// 1.e-6 - 1.e+14;  20*16 + 1 = 321
	nb_d = init_log_grid(delta_arr, d_min, d_max, point_nb_per_order);		// 1.e-2 - 1.e+9;   22*8 +1 = 177 

    ss.str("");
    ss << "./output/lvg_loss_func_q.txt";

    output.open(ss.str().c_str(), std::ios_base::out);
    output << scientific;
    output.precision(5);

    output << left << setw(18) << "";
    for (j = 0; j < nb_g; j++) {
        output << left << setw(13) << gamma_arr[j];
    }
    output << endl;

	// Calculation of the loss probability functions;
#pragma omp parallel private(i, j)
	{
		lvg_func_mu calc_q;
		calc_q.set_func(new lvg_func_cont_x1(1), new lvg_func_x2());

		double **p1 = alloc_2d_array<double>(nb_d, nb_g);
		memset(*p1, 0, nb_d*nb_g*sizeof(double));

#pragma omp for schedule(dynamic, 1), ordered
		for (i = 0; i < nb_d; i++) {
            for (j = 0; j < nb_g; j++)
		    {		
				calc_q.set_parameters(delta_arr[i], gamma_arr[j]);
				// one sided loss probability function for photon created in continuum processes is 0.5-0.5*q, where q is:
                p1[i][j] = qromb<lvg_func_mu>(calc_q, MU_MIN, 1., 1.e-5, false);
                
                // the photons with the mu < mu_min must be taken into account, ?
                p1[i][j] += MU_MIN*calc_q(MU_MIN);  
                p1[i][j] /= (SQRT_PI*delta_arr[i]);
			}

			#pragma omp ordered 
            {
				cout << delta_arr[i] << endl;    
                output << left << setw(5) << i << setw(13) << delta_arr[i];
                for (j = 0; j < nb_g; j++) {
                    output << left << setw(13) << p1[i][j];
                }
                output << endl;
			}
		}
		free_2d_array(p1);
	}
	output.close();

	delete [] gamma_arr;
	delete [] delta_arr;
}

void calc_line_overlap1(double d_min, double d_max, double g_min, double g_max, double gratio_min, double gratio_max, double dx_min, double dx_max, 
    int point_nb_per_order, int nb_dx_points, bool rect_profile)
{
    int i, j, k, l, m, nb_g, nb_gr, nb_d;
    double* gamma_arr, * gratio_arr, *dx_arr, *delta_arr;
    double** p_arr;

    // auxiliary arrays:
    int nb_g1, nb_g2, nb_g3;
    double* g1_arr, * g2_arr, * g3_arr;

    time_t timer = time(NULL);
    stringstream ss;
    ofstream output;

#ifdef __linux__
    ss.str("");
    ss << "./output/out_p1.txt";

    std::ofstream outerr;
    std::streambuf* orig_cerr;

    std::ofstream out;
    std::streambuf* orig_cout;

    outerr.open(ss.str().c_str(), std::ios::app);		// "/dev/null"
    orig_cerr = std::cerr.rdbuf(outerr.rdbuf());	// "/dev/null"

    out.open(ss.str().c_str(), std::ios::app);
    orig_cout = std::cout.rdbuf(out.rdbuf());
#endif

    cout << scientific;
    cout.precision(4);

    // Initialization of gamma array, 
    nb_g1 = init_log_grid(g1_arr, g_min, 1.e-2, 2);
    nb_g2 = init_log_grid(g2_arr, 1.e-2, 1.e-1, 4);
    nb_g3 = init_log_grid(g3_arr, 1.e-1, g_max, point_nb_per_order);

    nb_g = nb_g1 + nb_g2 + nb_g3 - 2;
    gamma_arr = new double [nb_g];
    
    memcpy(gamma_arr, g1_arr, nb_g1 * sizeof(double));
    memcpy(gamma_arr + nb_g1 - 1, g2_arr, nb_g2 * sizeof(double));
    memcpy(gamma_arr + nb_g1 + nb_g2 - 2, g3_arr, nb_g3 * sizeof(double));
  
    nb_d = init_log_grid(delta_arr, d_min, d_max, 2);   
    nb_gr = init_log_grid(gratio_arr, gratio_min, gratio_max, point_nb_per_order);

    init_lin_grid(dx_arr, dx_min, dx_max, nb_dx_points);
    p_arr = alloc_2d_array<double>(nb_gr, nb_g);

    ss.str("");
    ss << "./output/line_overlap_func_p1.txt";

    output.open(ss.str().c_str(), std::ios_base::out);
    output << scientific;
    output.precision(4);

    output << left << "# J1 = S1*K11 + S2*K12, here the K11 is tabulated," << endl;
    output << left << "# gamma (g1) is on horizontal axis, gamma ratio (g2/g1) - on vertical, for each delta and dx," << endl;
    output << left << "# number of points: delta (dust opacity), dx (frequency shift), gamma ratio (g2/g1), gamma of first line," << endl;
    output << left << setw(5) << nb_d << setw(5) << nb_dx_points << setw(5) << nb_gr << setw(5) << nb_g << endl;
    output << left << setw(17) << "";
    for (j = 0; j < nb_g; j++) {
        output << left << setw(12) << gamma_arr[j];
    }
    output << endl;

    // Calculation of the loss probability functions;
    for (l = 0; l < nb_d; l++) {
        for (k = 0; k < nb_dx_points; k++) {
            memset(*p_arr, 0, nb_gr * nb_g * sizeof(double));

#pragma omp parallel shared(l, k, nb_g, nb_gr, delta_arr, gamma_arr, gratio_arr, dx_arr) private(i, j, m)
            {
                line_overlap_rect1 calc_p_rect; // for rectangular profile;
                line_overlap_func_mu calc_p;
                
                line_profile_gauss* lprofile = new line_profile_gauss(X_RANGE_LINE_OVERLAP);
                
                line_overlap_sf1* f1 = new line_overlap_sf1();
                line_overlap_func2* f2 = new line_overlap_func2();

                f1->set_line_profile(lprofile);
                f2->set_line_profile(lprofile);
                f1->set_func(f2);
                
                calc_p.set_func(f1);
                calc_p.set_integr_lim(lprofile->x_min, lprofile->x_max);
                
                double** p1 = alloc_2d_array<double>(nb_gr, nb_g);
                memset(*p1, 0, nb_gr * nb_g * sizeof(double));

#pragma omp for schedule(dynamic, 1)
                for (m = 0; m < nb_g * nb_gr; m++) {
                    i = m / nb_g;
                    j = m % nb_g;
                    
                    // one sided loss probability function for photon created in line processes is 0.5-0.5*K, where K is p1
                    if (rect_profile) {
                        calc_p_rect.set_parameters(delta_arr[l], gamma_arr[j], gamma_arr[j] * gratio_arr[i], dx_arr[k]);
                        p1[i][j] = qromb<line_overlap_rect1>(calc_p_rect, MU_MIN, 1., 1.e-5, false) / gamma_arr[j];
                    }
                    else {
                        calc_p.set_parameters(delta_arr[l], gamma_arr[j], gamma_arr[j] * gratio_arr[i], dx_arr[k]);
                        p1[i][j] = qromb<line_overlap_func_mu>(calc_p, MU_MIN, 1., 1.e-5, false) / gamma_arr[j];
                    }              
                }
#pragma omp critical
                {
                    for (i = 0; i < nb_gr; i++) {
                        for (j = 0; j < nb_g; j++) {
                            p_arr[i][j] += p1[i][j];
                        }
                    }
                }
                free_2d_array(p1);
                delete lprofile;
                delete f1;
                delete f2;
            }
            output << left << setw(12) << delta_arr[l] << setw(12) << dx_arr[k] << endl;
            for (i = 0; i < nb_gr; i++) {
                output << left << setw(5) << i << setw(12) << gratio_arr[i];
                for (j = 0; j < nb_g; j++) {
                    output << left << setw(12) << p_arr[i][j];
                }
                output << endl;
            }
            cout << left << setw(12) << delta_arr[l] << setw(12) << dx_arr[k] << " time: " << (int)(time(NULL) - timer) << endl;
        }
    }
    output.close();

    delete[] g1_arr;
    delete[] g2_arr;
    delete[] g3_arr;

    delete[] delta_arr;
    delete[] gamma_arr;
    delete[] gratio_arr;
    delete[] dx_arr;
    free_2d_array(p_arr);
}

void calc_line_overlap2(double d_min, double d_max, double g_min, double g_max, double gratio_min, double gratio_max, double dx_min, double dx_max, 
    int point_nb_per_order, int nb_dx_points, bool rect_profile)
{
    int i, j, k, l, m, nb_g, nb_gr, nb_d;
    double* gamma_arr, * gratio_arr, * dx_arr, *delta_arr;
    double** p_arr;

    time_t timer = time(NULL);
    stringstream ss;
    ofstream output;

#ifdef __linux__
    ss.str("");
    ss << "./output/out_p2.txt";

    std::ofstream outerr;
    std::streambuf* orig_cerr;

    std::ofstream out;
    std::streambuf* orig_cout;

    outerr.open(ss.str().c_str(), std::ios::app);		// "/dev/null"
    orig_cerr = std::cerr.rdbuf(outerr.rdbuf());	// "/dev/null"

    out.open(ss.str().c_str(), std::ios::app);
    orig_cout = std::cout.rdbuf(out.rdbuf());
#endif

    cout << scientific;
    cout.precision(4);

    nb_g = init_log_grid(gamma_arr, g_min, g_max, point_nb_per_order);
    nb_d = init_log_grid(delta_arr, d_min, d_max, 2);
    nb_gr = init_log_grid(gratio_arr, gratio_min, gratio_max, point_nb_per_order);

    init_lin_grid(dx_arr, dx_min, dx_max, nb_dx_points);
    p_arr = alloc_2d_array<double>(nb_gr, nb_g);

    ss.str("");
    ss << "./output/line_overlap_func_p2.txt";

    output.open(ss.str().c_str(), std::ios_base::out);
    output << scientific;
    output.precision(4);

    output << left << "# J1 = S1*K11 + S2*K12, here the K12 is tabulated," << endl;
    output << left << "# gamma (g1) is on horizontal axis, gamma ratio (g2/g1) - on vertical, for each delta and dx," << endl;
    output << left << "# number of points: delta (dust opacity), dx (frequency shift), gamma ratio (g2/g1), gamma of first line," << endl;
    output << left << setw(5) << nb_d << setw(5) << nb_dx_points << setw(5) << nb_gr << setw(5) << nb_g << endl;
    output << left << setw(17) << "";
    for (j = 0; j < nb_g; j++) {
        output << left << setw(12) << gamma_arr[j];
    }
    output << endl;

    // Calculation of the loss probability functions;
    for (l = 0; l < nb_d; l++) {
        for (k = 0; k < nb_dx_points; k++) {
            memset(*p_arr, 0, nb_gr * nb_g * sizeof(double));

#pragma omp parallel shared(l, k, nb_g, nb_gr, delta_arr, gamma_arr, gratio_arr, dx_arr) private(i, j, m)
            {
                line_overlap_rect2 calc_p_rect; // for rectangular profile;
                line_overlap_func_mu calc_p;

                line_profile_gauss* lprofile = new line_profile_gauss(X_RANGE_LINE_OVERLAP);
                
                line_overlap_sf2* f1 = new line_overlap_sf2();
                line_overlap_func2* f2 = new line_overlap_func2();

                f1->set_line_profile(lprofile);
                f2->set_line_profile(lprofile);
                f1->set_func(f2);

                calc_p.set_func(f1);
                calc_p.set_integr_lim(lprofile->x_min + dx_arr[k], lprofile->x_max + dx_arr[k]);

                double** p1 = alloc_2d_array<double>(nb_gr, nb_g);
                memset(*p1, 0, nb_gr * nb_g * sizeof(double));

#pragma omp for schedule(dynamic, 1)
                for (m = 0; m < nb_g * nb_gr; m++) {
                    i = m / nb_g;
                    j = m % nb_g;

                    // Note!!! The correction was made - 1/g2 instead of 1./g1
                    // one sided loss probability function for photon created in line processes is 0.5-0.5*K, where K is 
                    if (rect_profile) {
                        calc_p_rect.set_parameters(delta_arr[l], gamma_arr[j], gamma_arr[j] * gratio_arr[i], dx_arr[k]);
                        p1[i][j] = qromb<line_overlap_rect2>(calc_p_rect, MU_MIN, 1., 1.e-5, false) / (gamma_arr[j] * gratio_arr[i]);
                    }
                    else {
                        calc_p.set_parameters(delta_arr[l], gamma_arr[j], gamma_arr[j] * gratio_arr[i], dx_arr[k]);
                        p1[i][j] = qromb<line_overlap_func_mu>(calc_p, MU_MIN, 1., 1.e-5, false) / (gamma_arr[j] * gratio_arr[i]);
                    }
                }
#pragma omp critical
                {
                    for (i = 0; i < nb_gr; i++) {
                        for (j = 0; j < nb_g; j++) {
                            p_arr[i][j] += p1[i][j];
                        }
                    }
                }
                delete f1;
                delete f2;
                delete lprofile;
                free_2d_array(p1);
            }
            output << left << setw(12) << delta_arr[l] << setw(12) << dx_arr[k] << endl;
            for (i = 0; i < nb_gr; i++) {
                output << left << setw(5) << i << setw(12) << gratio_arr[i];
                for (j = 0; j < nb_g; j++) {
                    output << left << setw(12) << p_arr[i][j];
                }
                output << endl;
            }
            cout << left << setw(12) << delta_arr[l] << setw(12) << dx_arr[k] << " time: " << (int)(time(NULL) - timer) << endl;
        }
    }
    output.close();

    delete[] delta_arr;
    delete[] gamma_arr;
    delete[] gratio_arr;
    delete[] dx_arr;
    free_2d_array(p_arr);
}

/*
void process_initialisation(int &ProcessNb);
void process_initialisation(int &ProcessNb)
{
	int l, nb, seconds; 
	time_t c_end;
	std::fstream iofile;
	
	srand(time(0) + getpid());
	seconds = rand()%MAX_TIME_DELAY;
	c_end = seconds + time(NULL);
		
	while (c_end > time(NULL));

	iofile.open("process_nb.txt", std::ios_base::in | std::ios_base::out);
	if (iofile.is_open())
	{
		iofile.seekg(0, std::ios::beg);
		iofile >> nb;
		iofile >> ProcessNb;

		iofile.seekp(0, std::ios::beg);
		iofile << nb << endl;
		for (l = ProcessNb+1; l <= nb; l++) {
			iofile << l << endl;
		}
		iofile.close();
	}
	else ProcessNb = 0;
}

// Derivative on delta:
class lvg_func_x2_dd : public lvg_func_x2
{
public:
	double operator() (double x2) {
		return exp(-x2*x2) *exp(-0.5*(ef.f(x2)-y1)/gz - (x2-x1)/dz) *(x2-x1)/dz; // 1/(d sqrt(pi)) is outside the integral;
	};
	bool condition(double t) { 
		return (exp(-0.5*(ef.f(x1+t)-y1)/gz - t/dz) < DBL_EPSILON); 
	};
	double approximation() {
		e = 2.*x1*gz + ONEDIVBY_SQRT_PI*exp(-x1*x1) + g/d;
		return gz*g*exp(-x1*x1)/(e*e*d);	// 1/(d sqrt(pi)) is outside the integral;
	};
};

// Derivative on g:
class lvg_func_x2_dg : public lvg_func_x2
{
public:
	double operator() (double x2) {
		e = 0.5*(ef.f(x2)-y1)/gz;
		return exp(-x2*x2) *exp(-e-(x2-x1)/dz) *(e-1.);	// 1/(g sqrt(pi)) is outside the integral;
	};
	bool condition(double t) { 
		return (exp(-0.5*(ef.f(x1+t)-y1)/gz - t/dz ) < DBL_EPSILON); 
	};
	double approximation() {
		e = 2.*x1*dz + ONEDIVBY_SQRT_PI*exp(-x1*x1) + g/d;
		return -exp(-x1*x1)*gz *(2.*gz*x1 + g/d) /(e*e); // 1/(g sqrt(pi)) is outside the integral;
	};
};

// Derivative on g and d
class lvg_func_x2_ddg : public lvg_func_x2
{
public:
	double ee;
	double operator() (double x2) {
		e = 0.5*(ef.f(x2)-y1)/gz;
		ee = (x2-x1)/dz;
		return exp(-x2*x2) *exp(-e-ee) *(e-1.) *ee; // 1/(g d sqrt(pi)) is outside the integral;
	};
	bool condition(double t) { 
		return (exp(-0.5*(ef.f(x1+t)-y1)/gz - t/dz) < DBL_EPSILON); 
	};
	double approximation() {
		e = 2.*x1*gz + ONEDIVBY_SQRT_PI*exp(-x1*x1) + g/d;
		return exp(-x1*x1)*gz*gz *(2.*ONEDIVBY_SQRT_PI*exp(-x1*x1)-e) /(dz*e*e*e); // 1/(g d sqrt(pi)) is outside the integral;
	};
};

void test_integration()
{
	double delta, gamma;
	lvg_func_mu calc_p, calc_p_dd, calc_p_dg, calc_p_ddg;

	calc_p.set_func(new lvg_func_x2());
	calc_p_dd.set_func(new lvg_func_x2_dd());
	calc_p_dg.set_func(new lvg_func_x2_dg());
	calc_p_ddg.set_func(new lvg_func_x2_ddg());

	delta = 1.0000e-02; 
	gamma = 7.4989e+09;
	
	calc_p_dd.set_parameters(delta, gamma);
	cout << qromb<lvg_func_mu>(calc_p_dd, 1.e-8, 1., 1.e-5) /(M_PI*gamma) << endl;
	getchar();
}

void calc_loss_func_deriv()
{
	int i, j, k, n, start, nb_g, nb_d;
	double min_gamma, max_gamma, min_delta, max_delta;
	double *gamma_arr, *delta_arr;
	double **p_dd_arr, **p_dg_arr, **p_ddg_arr;
	
	stringstream ss;
	ofstream output;

	lvg_func_mu calc_p_dd, calc_p_dg, calc_p_ddg;

	cout << scientific;
	cout.precision(6);

	//k = 11;
	process_initialisation(k);
	n = 2;
	start = n*k + 17*4 + 8*2;

#ifdef __linux__
	ss.str("");
	ss << "./output/out_d";
	ss << k;
	ss << ".txt";

	std::ofstream outerr;
	std::streambuf *orig_cerr;
	
	std::ofstream out;
	std::streambuf *orig_cout; 

	outerr.open(ss.str().c_str(), std::ios::app);		// "/dev/null"
	orig_cerr = std::cerr.rdbuf(outerr.rdbuf());	// "/dev/null"
	
	out.open(ss.str().c_str(), std::ios::app);
	orig_cout = std::cout.rdbuf(out.rdbuf());
#endif
	
	nb_g = init_log_grid(gamma_arr, min_gamma = 1.e-6, max_gamma = 1.e+14, 16);		// 20*16 + 1 = 321
	nb_d = init_log_grid(delta_arr, min_delta = 1.e-2, max_delta = 1.e+9, 16);		// 22*8 + 1 = 177

	calc_p_dd.set_func(new lvg_func_x2_dd());
	calc_p_dg.set_func(new lvg_func_x2_dg());
	calc_p_ddg.set_func(new lvg_func_x2_ddg());

	p_dd_arr = alloc_2d_array<double>(nb_d, nb_g);
	memset(*p_dd_arr, 0, nb_d*nb_g*sizeof(double));

	p_dg_arr = alloc_2d_array<double>(nb_d, nb_g);
	memset(*p_dg_arr, 0, nb_d*nb_g*sizeof(double));

	p_ddg_arr = alloc_2d_array<double>(nb_d, nb_g);
	memset(*p_ddg_arr, 0, nb_d*nb_g*sizeof(double));

	// Calculation of the loss probability functions;
	for (i = start; (i < start+n) && (i < nb_d); i++)
	{	
		for (j = 0; j < nb_g; j++)
		{	
			calc_p_dd.set_parameters(delta_arr[i], gamma_arr[j]);
			p_dd_arr[i][j] = qromb<lvg_func_mu>(calc_p_dd, 1.e-8, 1., 1.e-5) /(M_PI*gamma_arr[j]*delta_arr[i]);
			
			cout << left << setw(5) << i << setw(5) << j << setw(15) << p_dd_arr[i][j];
			
			calc_p_dg.set_parameters(delta_arr[i], gamma_arr[j]);
			p_dg_arr[i][j] = qromb<lvg_func_mu>(calc_p_dg, 1.e-8, 1., 1.e-5) /(M_PI*gamma_arr[j]*gamma_arr[j]);
			
			cout << left << setw(15) << p_dg_arr[i][j];
				
			calc_p_ddg.set_parameters(delta_arr[i], gamma_arr[j]);
			p_ddg_arr[i][j] = qromb<lvg_func_mu>(calc_p_ddg, 1.e-8, 1., 1.e-5) /(M_PI*gamma_arr[j]*gamma_arr[j]*delta_arr[i]);

			cout << left << setw(15) << p_ddg_arr[i][j] << endl;
		}
	}

	ss.str("");
	ss << "./output/lvg_loss_func_deriv";
	ss << k;
	ss << ".txt";

	output.open(ss.str().c_str(), std::ios_base::out);
	output << scientific;
	output.precision(6);

	for (j = 0; j < nb_g; j++) {
		output << left << setw(15) << gamma_arr[j];
	}
	output << endl;

	for (i = start; (i < start+n) && (i < nb_d); i++) 
	{
		output << left << setw(5) << i << setw(15) << delta_arr[i];
		for (j = 0; j < nb_g; j++) {
			output << left << setw(15) << p_dd_arr[i][j];
		}
		output << endl;
	}
	output << endl;

	for (i = start; (i < start+n) && (i < nb_d); i++) 
	{
		output << left << setw(5) << i << setw(15) << delta_arr[i];
		for (j = 0; j < nb_g; j++) {
			output << left << setw(15) << p_dg_arr[i][j];
		}
		output << endl;
	}
	output << endl;

	for (i = start; (i < start+n) && (i < nb_d); i++) 
	{
		output << left << setw(5) << i << setw(15) << delta_arr[i];
		for (j = 0; j < nb_g; j++) {
			output << left << setw(15) << p_ddg_arr[i][j];
		}
		output << endl;
	}
	output.close();

	delete [] gamma_arr;
	delete [] delta_arr;

	free_2d_array(p_dd_arr);
	free_2d_array(p_dg_arr);
	free_2d_array(p_ddg_arr);
}

void test_func_deriv()
{
	const int text_width = 160;
	char text_line[text_width];
	int i, j, nb_d, nb_g;
	double deriv;
	double *delta_arr, *gamma_arr;
	double **p_arr, **p_d1_arr, **p_d2_arr, **p_d12_arr;
	
	string file_name;
	stringstream ss;
	ifstream input;
	ofstream output;

	file_name = "lvg_loss_func.txt";
	input.open(file_name.c_str(), ios_base::in);

	if (!input.is_open()) {
		cout << "Can't open " << file_name << endl;
		exit(1);
	}

	input.getline(text_line, text_width-1);
	input.getline(text_line, text_width-1);
	input.getline(text_line, text_width-1);

	input >> nb_d >> nb_g;

	delta_arr = new double [nb_d];
	gamma_arr = new double [nb_g];
	
	p_arr = alloc_2d_array<double>(nb_d, nb_g);
	p_d1_arr = alloc_2d_array<double>(nb_d, nb_g);
	p_d2_arr = alloc_2d_array<double>(nb_d, nb_g);
	p_d12_arr = alloc_2d_array<double>(nb_d, nb_g);

	for (i = 0; i < nb_g; i++) {
		input >> gamma_arr[i];
	}

	for (i = 0; i < nb_d; i++) 
	{
		input >> j >> delta_arr[i];
		for (j = 0; j < nb_g; j++) {
			input >> p_arr[i][j];
		}
	}
	input.close();

// The estimation of the derivatives; derivative on delta:
	for (j = 0; j < nb_g; j++) 
	{
		for (i = 1; i < nb_d-1; i++) 
		{
			p_d1_arr[i][j] = 0.5*( (p_arr[i+1][j] - p_arr[i][j])/(delta_arr[i+1] - delta_arr[i]) +
				(p_arr[i][j] - p_arr[i-1][j])/(delta_arr[i] - delta_arr[i-1]) );		
		}
		p_d1_arr[0][j] = (p_arr[1][j] - p_arr[0][j])/(delta_arr[1] - delta_arr[0]);
		p_d1_arr[nb_d-1][j] = (p_arr[nb_d-1][j] - p_arr[nb_d-2][j])/(delta_arr[nb_d-1] - delta_arr[nb_d-2]);
	}		

	// derivative on gamma:
	for (i = 0; i < nb_d; i++) 
	{
		for (j = 1; j < nb_g-1; j++) 
		{
			p_d2_arr[i][j] = 0.5*( (p_arr[i][j+1] - p_arr[i][j])/(gamma_arr[j+1] - gamma_arr[j]) +
				(p_arr[i][j] - p_arr[i][j-1])/(gamma_arr[j] - gamma_arr[j-1]) ); 
		}
		p_d2_arr[i][0] = (p_arr[i][1] - p_arr[i][0])/(gamma_arr[1] - gamma_arr[0]);
		p_d2_arr[i][nb_g-1] = (p_arr[i][nb_g-1] - p_arr[i][nb_g-2])/(gamma_arr[nb_g-1] - gamma_arr[nb_g-2]);
	}
	
	// cross-derivative:
	for (i = 0; i < nb_d; i++) 
	{
		for (j = 1; j < nb_g-1; j++) 
		{
			p_d12_arr[i][j] = 0.25*( (p_d1_arr[i][j+1] - p_d1_arr[i][j])/(gamma_arr[j+1] - gamma_arr[j]) +
				(p_d1_arr[i][j] - p_d1_arr[i][j-1])/(gamma_arr[j] - gamma_arr[j-1]) ); 
		}
		p_d12_arr[i][0] = 0.5*(p_d1_arr[i][1] - p_d1_arr[i][0])/(gamma_arr[1] - gamma_arr[0]);
		p_d12_arr[i][nb_g-1] = 0.5*(p_d1_arr[i][nb_g-1] - p_d1_arr[i][nb_g-2])/(gamma_arr[nb_g-1] - gamma_arr[nb_g-2]);
	}

	for (j = 0; j < nb_g; j++) 
	{
		for (i = 1; i < nb_d-1; i++) 
		{
			p_d12_arr[i][j] += 0.25*( (p_d2_arr[i+1][j] - p_d2_arr[i][j])/(delta_arr[i+1] - delta_arr[i]) +
				(p_d2_arr[i][j] - p_d2_arr[i-1][j])/(delta_arr[i] - delta_arr[i-1]) ); 
		}
		p_d12_arr[0][j] += 0.5*(p_d2_arr[1][j] - p_d2_arr[0][j])/(delta_arr[1] - delta_arr[0]);
		p_d12_arr[nb_d-1][j] += 0.5*(p_d2_arr[nb_d-1][j] - p_d2_arr[nb_d-2][j])/(delta_arr[nb_d-1] - delta_arr[nb_d-2]);
	}

	file_name = "lvg_loss_func_deriv.txt";
	input.open(file_name.c_str(), ios_base::in);

	if (!input.is_open()) {
		cout << "Can't open " << file_name << endl;
		exit(1);
	}

	input.getline(text_line, text_width-1);
	input.getline(text_line, text_width-1);
	input.getline(text_line, text_width-1);

	input >> nb_d >> nb_g;
	
	for (i = 0; i < nb_g; i++) {
		input >> gamma_arr[i];
	}

	for (i = 0; i < nb_d; i++) 
	{
		input >> j >> delta_arr[i];
		for (j = 0; j < nb_g; j++) {
			input >> deriv;
			p_d1_arr[i][j] = fabs((p_d1_arr[i][j] - deriv)/(deriv + DBL_EPSILON));
		}
	}

	for (i = 0; i < nb_d; i++) 
	{
		input >> j >> delta_arr[i];
		for (j = 0; j < nb_g; j++) {
			input >> deriv;
			p_d2_arr[i][j] = fabs((p_d2_arr[i][j] - deriv)/(deriv + DBL_EPSILON));
		}
	}

	for (i = 0; i < nb_d; i++) 
	{
		input >> j >> delta_arr[i];
		for (j = 0; j < nb_g; j++) {
			input >> deriv;
			p_d12_arr[i][j] = fabs((p_d12_arr[i][j] - deriv)/(deriv + DBL_EPSILON));
		}
	}
	input.close();

	ss.str("");
	ss << "./output/lvg_func_deriv_comp.txt";

	output.open(ss.str().c_str(), std::ios_base::out);
	output << scientific;
	output.precision(6);

	for (i = 0; i < nb_d; i++) {
		for (j = 0; j < nb_g; j++) {
			output << left << setw(5) << i << setw(5) << j << setw(15) << p_d1_arr[i][j] 
				<< setw(15) << p_d2_arr[i][j] << setw(15) << p_d12_arr[i][j] << endl;
		}
	}
	output.close();

	delete [] gamma_arr;
	delete [] delta_arr;

	free_2d_array(p_arr);
	free_2d_array(p_d1_arr);
	free_2d_array(p_d2_arr);
	free_2d_array(p_d12_arr);
}
*/
