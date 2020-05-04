#ifndef GAMMA_JUMP
#define GAMMA_JUMP

#include <stan/math/prim/mat.hpp>
#include "base_jump.hpp"
#include "../rng.hpp"

using namespace stan::math;

class GammaJump: public BaseJump {
 protected:
    double alpha, beta;

  public:
    GammaJump(double alpha, double beta): alpha(alpha), beta(beta) {}
    ~GammaJump() {}
    
    double sample_tilted(double u) override;

    double sample_given_data(
        int ndata, double curr, double u) override;

    double laplace(double u) override;
};

#endif