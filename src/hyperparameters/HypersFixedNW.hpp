#ifndef HYPERSFIXEDNW_HPP
#define HYPERSFIXEDNW_HPP

#include <cassert>
#include <Eigen/Dense>

class HypersFixedNW {
private:
    Eigen::VectorXd mu0;
    Eigen::MatrixXd lambda0;

public:
    // Destructor and constructor
    ~HypersFixedNW() = default;

    HypersFixedNW(const Eigen::VectorXd &mu0_, const Eigen::MatrixXd &lambda0_):
        mu0(mu0_), lambda0(lambda0_){}

    // Getters
    const Eigen::VectorXd get_mu0(){return mu0;}
    const Eigen::MatrixXd get_lambda0(){return lambda0;}

    // Setters
    void set_mu0(const Eigen::VectorXd &mu_0_){mu0 = mu_0_;}
    void set_lambda0(const Eigen::MatrixXd &lambda_0_){lambda0 = lambda_0_;}
};


#endif // HYPERSFIXEDNW_HPP
