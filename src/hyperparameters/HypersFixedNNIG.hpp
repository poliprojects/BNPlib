#ifndef HYPERSFIXEDNNIG_HPP
#define HYPERSFIXEDNNIG_HPP

#include <cassert>

class HypersFixedNNIG {
private:
    double mu0, lambda, alpha0, beta0;

public:
    // Destructor and constructor
    ~HypersFixedNNIG() = default;

    HypersFixedNNIG(const double mu0_, const double lambda_,
        const double alpha0_, const double beta0_):
        mu0(mu0_), lambda(lambda_), alpha0(alpha0_), beta0(beta0_) {
        assert(lambda > 0);
        assert(alpha0 > 0);
        assert(beta0  > 0);
    }

    // Getters
    const double get_mu0(){return mu0;}
    const double get_alpha0(){return alpha0;}
    const double get_beta0(){return beta0;}
    const double get_lambda(){return lambda;}

    // Setters
    void set_mu0(const double mu_0_){mu0 = mu_0_;}
    void set_alpha0(const double alpha_0_){alpha0 = alpha_0_;}
    void set_beta0(const double beta_0_){beta0 = beta_0_;}
    void set_lambda(const double lambda_){lambda = lambda_;}
};


#endif // HYPERSFIXEDNNIG_HPP
