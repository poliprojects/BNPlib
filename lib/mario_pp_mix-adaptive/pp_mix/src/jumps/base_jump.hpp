#ifndef BASE_JUMP_HPP
#define BASE_JUMP_HPP

#include <vector>
#include <Eigen/Dense>

using namespace Eigen;

class BaseJump {
 public:
    virtual ~BaseJump() {};

    virtual double sample_tilted(double u) = 0;

    virtual double sample_given_data(
        int ndata, double curr, double u) = 0;

    virtual double laplace(double u) = 0;
};


#endif 