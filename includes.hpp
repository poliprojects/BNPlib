#ifndef INCLUDES_HPP
#define INCLUDES_HPP

// Universal includes, to be included in every single header:
#include <Eigen/Dense> 
#include <vector>
#include <armadillo>
#include <stan/math/prim/mat.hpp>
#include <boost/random/random_number_generator.hpp>
#include <boost/random/detail/qrng_base.hpp>

namespace bnplib{
    using data_t = double;
    using par_t = double;
    using parvec_t = std::vector<par_t>;
}

using namespace bnplib;

#endif // INCLUDES_HPP
