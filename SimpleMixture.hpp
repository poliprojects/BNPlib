#ifndef SIMPLEMIXTURE_HPP
#define SIMPLEMIXTURE_HPP

#include "includes_universal.hpp"
//#include <Eigen/Dense>
//#include <stan/math/prim/mat.hpp>
//#include <boost/random/random_number_generator.hpp>
//#include <boost/random/detail/qrng_base.hpp>

class SimpleMixture {

private:
    double totalmass;

public:
    ~SimpleMixture() = default;

    SimpleMixture(double totalmass): totalmass(totalmass) {
      assert(totalmass>=0);
    }

    double const get_totalmass(){return totalmass;}

};

#endif // SIMPLEMIXTURE_HPP
