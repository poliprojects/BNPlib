#ifndef HYPERSDUMMY_HPP
#define HYPERSDUMMY_HPP

#include <cassert>
#include <Eigen/Dense>

class HypersDummy {
private:
    Eigen::VectorXd mu0;
    Eigen::MatrixXd lambda0;

public:
    // Destructor and constructor
    ~HypersDummy() = default;

    HypersDummy(const Eigen::VectorXd &mu0, const Eigen::MatrixXd &lambda0):
        mu0(mu0), lambda0(lambda0){}

    // Getters
    const Eigen::VectorXd get_mu0(){return mu0;}
    const Eigen::MatrixXd get_lambda0(){return lambda0;}

    // Setters
    void set_mu0(const Eigen::VectorXd &mu0){mu0 = mu0;}
    void set_lambda0(const Eigen::MatrixXd &lambda0){lambda0 = lambda0;}
};


#endif // HYPERSDUMMY_HPP
