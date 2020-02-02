#ifndef INCLUDES_UNIVERSAL_HPP
#define INCLUDES_UNIVERSAL_HPP

#include <Eigen/Dense> 
// typedef Matrix<double, Dynamic, 1> VectorXd

//#include <random>

namespace bnplib{
    using data_t = double;
    using par_t = double;
    using parvec_t = std::vector<par_t>;
}

using namespace bnplib;

#endif // INCLUDES_UNIVERSAL_HPP
