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

    // Setters
    void set_mu0(const Eigen::VectorXd &mu0_){
    	assert(mu0_.size() == mu0.size());
    	mu0 = mu0_;
    }
    void set_lambda(const double lambda_){
    	assert(lambda_ > 0);
    	lambda = lambda_;
    }
    void set_tau0(const Eigen::MatrixXd &tau0_) {
    	assert(tau0_.rows() == tau0_.cols());
    	assert(mu0.size() == tau0_.rows());
    	tau0 = tau0_;
    }
    void set_nu(const double nu_){
    	assert(nu_ > mu0.size()-1);
    	nu = nu_;
    }
};


#endif // HYPERSFIXEDNW_HPP
