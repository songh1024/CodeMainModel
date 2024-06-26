#include <iostream>
#include <string>
#include <iomanip>
#include <cstdio>
#include <fstream>
#include <math.h>
#include "function_declaration.h" 

#ifndef PARAMETERS_H_   // To make sure I don't declare the function more than once by including the header multiple times.
#define PARAMETERS_H_

const int lag = 10;

const double beta = 0.89;
const double sigma = 1.5;

const double shape = 4.7;
const double scale = 1.0;
const double cutoff = 0.99;
const double grammar = 0.28;

const double time_preference_rate = 1 / beta - 1;

const double delta = 0.08;
const double alpha = 0.33;
const double nu = 0.14;

const double ksi = 0.59;
const double chi = 0.049;
const double conspsi = 0.83;
const double thetapsi = 0.215;
const double conszeta = 0.1;
const double thetazeta = 0.0015;

const double kappa_ini = 1e10;
const double eta_ini = 0.0;
const double kappa_new = 2.2;
const double eta_new = 0.008;

const int nb = 201;
const int ntheta = 6;
const int nz = 7;

const double minb = 0.0;
const double maxb = 200.0;
const double mintheta = 0.0;
const double maxtheta = 1.0;

const int maxit_DP = 100;
const int maxit_gss = 200;
const double tol_itvfun = 1e-4;
const double tol_gss = 1e-4;
const double golden_ratio = (pow(5, 0.5) - 1) / 2;

const int coefficient = ntheta;

const int writing_precision = 6;

const int nb_pdf = 10001;
const int ntheta_pdf = 21;

const double initial_saving = 2.0;
const int period = 201;

const int nmkt = 1428;
const double f_w = 0.769;

const double h = 150;

const int nyears = 11;
const int period_trans = 200;

const int nincome_pdf = 100001;
const double min_income = 0.0;
const double max_income = 50.0;

extern double kappa;
extern double eta;

#endif