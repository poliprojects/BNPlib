#ifndef INCLUDES_UNIVERSAL_HPP
#define INCLUDES_UNIVERSAL_HPP

// Universal includes, to be included in every single header:
#include <vector>
#include <armadillo>
#include <assert.h> // assert
#include <memory> // smart pointers
//#include <Eigen/Dense> // TODO why does it not work here?
//#include <stan/math/prim/mat.hpp> // TODO same
//#include <boost/random/random_number_generator.hpp> // TODO same
//#include <boost/random/detail/qrng_base.hpp> // TODO same

namespace bnplib{
    using data_t = double;
    using par_t = double;
    using parvec_t = std::vector<par_t>;
}

using namespace bnplib;

#endif // INCLUDES_UNIVERSAL_HPP
