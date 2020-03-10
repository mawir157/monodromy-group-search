#pragma once

#include <algorithm>
#include <armadillo>
#include <complex>
#include <vector>
#include <iostream> 
#include <set>
#include <string>
#include <cmath>
#include <chrono>
#include <ctime>
#include <fstream>
#include <random>

#include <math.h> 


typedef std::complex<double> comp_d;
typedef arma::Mat<comp_d> CompMat3;
typedef arma::cx_dvec Point;
typedef arma::cx_dvec Polar;

static const double PI = 3.14159265358979323;
static const comp_d omega(-1.0 / 2.0, std::sqrt(3.0) / 2.0);
static const double LOWER_TOL = 1e-10;
static const double TOL = 1e-6;
static const unsigned int MAX_ORDER = 1000;
static const unsigned int MAX_BRAID = 1000;
static const unsigned int MAX_WORD  = 65536;
