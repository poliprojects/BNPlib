#ifndef HYPERSFIXEDNW_HPP
#define HYPERSFIXEDNW_HPP

#include <cassert>
#include <Eigen/Dense>

class HypersFixedNW {
private:
    Eigen::VectorXd mu0;
    double lambda;
    Eigen::MatrixXd tau0;
    double nu;

    Eigen::MatrixXd tau0_inv; // TODO etc: serve?

public:
    // Destructor and constructor
    ~HypersFixedNW() = default;

    HypersFixedNW(const Eigen::VectorXd &mu0_, const double lambda_,
        const Eigen::MatrixXd &tau0_, const double nu_):
        mu0(mu0_), lambda(lambda_), tau0(tau0_), nu(nu_) {
            // Check validity of parameters
            unsigned int dim = mu0.size();
            assert(lambda > 0);
            assert(dim == tau0.rows());
            assert(tau0.rows() == tau0.cols());
            assert(nu > dim-1);
            // TODO assert tau0 pos def matrix?

            tau0_inv = tau0.inverse();
        }

    // Getters
    const Eigen::VectorXd get_mu0(){return mu0;}
    const double get_lambda(){return lambda;}
    const Eigen::MatrixXd get_tau0(){return tau0;}
    const double get_nu(){return nu;}
    const Eigen::MatrixXd get_tau0_inv(){return tau0_inv;}

    // Setters
    void set_mu0(const Eigen::VectorXd &mu0_){mu0 = mu0_;}
    void set_lambda(const double lambda_){lambda = lambda_;}
    void set_tau0(const Eigen::MatrixXd &tau0_) {
        tau0 = tau0_;
        tau0_inv = tau0.inverse();
    }
    void set_nu(const double nu_){nu = nu_;}
    // TODO checks like in the constructor
};


#endif // HYPERSFIXEDNW_HPP
