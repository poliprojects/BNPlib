//
// Created by mario on 30/04/19.
//

#ifndef CPPMODEL_OPTIONS_HPP
#define CPPMODEL_OPTIONS_HPP

#include <assert.h>
#include <algorithm>
#include <iostream>
#include <limits>
#include <map>
#include <math.h>
#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#ifdef USE_RCPP_ARMADILLO
  #include <RcppArmadillo.h>
#else
  #include <armadillo>
#endif

#include "statslib/stats.hpp"

using std::vector;

static double infinity = std::numeric_limits<float>::infinity();

static double EPS = 1e-10;


#endif //CPPMODEL_OPTIONS_HPP
