#ifndef HYPERSFIXED_HPP
#define HYPERSFIXED_HPP

class HypersFixed {

private:
    double mu0, lambda, alpha0, beta0;

public:
    double get_mu0(){return mu0;}
    double get_alpha0(){return alpha0;}
    double get_beta0(){return beta0;}
    double get_lambda(){return lambda;}

    ~HypersFixed() = default;

    HypersFixed(double mu0, double lambda, double alpha0, double beta0):
        mu0(mu0), lambda(lambda), alpha0(alpha0), beta0(beta0) {
            assert(lambda > 0);
            assert(alpha0 > 0);
            assert(beta0  > 0);
        }
};

#endif // HYPERSFIXED_HPP
