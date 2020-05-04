#ifndef SIMULATE_STRAUSS_HPP
#define SIMULATE_STRAUSS_HPP

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "pointprocess/point2pattern.h"
#include "pointprocess/pointprocess.h"
#include "pointprocess/sampler.h"
#include "pointprocess/strauss.h"

#include <stdlib.h>
#include <Eigen/Dense>

Eigen::MatrixXd simulate_strauss_moller(
    const Eigen::MatrixXd &ranges, double beta, double gamma, double R);

#endif 