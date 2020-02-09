#ifndef HYPERSFIXEDNNIG_HPP
#define HYPERSFIXEDNNIG_HPP

class HypersFixedNNIG {

private:
    double mu0, lambda, alpha0, beta0;

public:
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
    void set_mu0(const double mu0_){mu0 = mu0_;}
    void set_alpha0(const double alpha0_){alpha0 = alpha0_;}
    void set_beta0(const double beta0_){beta0 = beta0_;}
    void set_lambda(const double lambda_){lambda = lambda_;}
};

#endif // HYPERSFIXEDNNIG_HPP
