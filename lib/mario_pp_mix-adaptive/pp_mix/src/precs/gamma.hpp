#ifndef PREC_GAMMA_HPP
#define PREC_GAMMA_HPP

#include "base_prec.hpp"
#include "../rng.hpp"
#include "../utils.hpp"

#include <algorithm>
#include <stan/math/prim/mat.hpp>

class GammaPrec: public BaseUnivPrec {
protected:
    double alpha;
    double beta;

public:
    GammaPrec(double alpha, double beta);

    ~GammaPrec() {}

    double sample_prior() override;

    double sample_given_data(
        const std::vector<double> &data, const double &curr,
        const VectorXd &mean) override;
};


#endif