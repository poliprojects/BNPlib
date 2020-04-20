#ifndef HYPERSFIXEDNNIG_HPP
#define HYPERSFIXEDNNIG_HPP

#include <cassert>

class HypersFixedNNIG {
private:
    double mu0, lambda, alpha0, beta0;

public:
    // Destructtor and constructor
    ~HypersFixedNNIG() = default;

    HypersFixedNNIG(const double mu0, const double lambda, const double alpha0,
        const double beta0):
        mu0(mu0), lambda(lambda), alpha0(alpha0), beta0(beta0) {
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
    void set_mu0(const double mu0){mu0 = mu0;}
    void set_alpha0(const double alpha0){alpha0 = alpha0;}
    void set_beta0(const double beta0){beta0 = beta0;}
    void set_lambda(const double lambda){lambda = lambda;}
};


#endif // HYPERSFIXEDNNIG_HPP
