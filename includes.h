#pragma once

#include <algorithm>
#include <armadillo>
#include <chrono>
#include <cmath>
#include <complex>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include <math.h>

typedef std::complex<double> comp_d;
typedef arma::Mat<comp_d> CompMat3;
typedef std::shared_ptr<const CompMat3> p_CompMat3;
typedef arma::cx_dvec Point;
typedef arma::cx_dvec Polar;

static const double PI = 3.14159265358979323;
static const comp_d omega(-1.0 / 2.0, std::sqrt(3.0) / 2.0);
static const comp_d c_zero(0.0, 0.0);
static const double LOWER_TOL = 1e-10;
static const double TOL = 1e-6;
static const double HIGHER_TOL = 1e-3;
static const unsigned int MAX_ORDER = 150;
static const unsigned int MAX_BRAID = 50;

enum RunMode { file, loop, verbose, matrix, unknown };
