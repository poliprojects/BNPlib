#ifndef NNIGHIERARCHY_HPP
#define NNIGHIERARCHY_HPP

#include "includes_universal.hpp"

// Normal likelihoood, Normal Inverse Gamma hierarchy
template<class Hypers> //Hypers = TupleWrapper, distro, ...
class NNIGHierarchy {
protected:
    using state_t = std::array<par_t,2>;

    std::mt19937 rng;
    state_t state; // current values for F's parameters: mu, sigma
    std::shared_ptr<Hypers> hypers; // current values for G0's parameters:
                                    // mu_0, lambda0, alpha, beta

public:
    // Constructors:
    ~NNIGHierarchy() = default;
    NNIGHierarchy(std::shared_ptr<Hypers> hypers): hypers(hypers) {}

    // Getters/setters:
    state_t get_state(){return state;}
    void set_state(const state_t &s){state = s;}
    void set_state(int pos, par_t val){state[pos] = val;}

    int get_count(){return hypers.use_count();}

    double log_like(data_t datum);

    void draw();

    void sample_given_data(std::vector<data_t> data);



    parvec_t normal_gamma_update(std::vector<data_t> data, double mu0,
        double alpha0, double beta0, double lambda0) {
        
        double mu_post, alpha_post, beta_post, lambda_post;
        int n = data.size();
        if (n == 0)
            return parvec_t{mu0, alpha0, beta0, lambda0};

        // Compute sample mean
        double y_bar = accumulate(data.begin(), data.end(), 0.0) / n;

        mu_post = (lambda0 * mu0 + n * y_bar) / (lambda0 + n);
        alpha_post = 1.0 * alpha0 + 1.0 * n / 2;

        // double ss = n * arma::var(data, 1); // divides by n, not n-1
        double ss = std::inner_product( data.begin(), data.end(), data.begin(),
            0.0 ) - n*y_bar*y_bar; // TODO check

        beta_post = (beta0 + 0.5 * ss + 0.5 * lambda0 /
            (n + lambda0) * n * std::pow((y_bar - mu0), 2));

        lambda_post = lambda0 + n;

        return parvec_t{mu_post, alpha_post, beta_post, lambda_post};
    }

};

# endif // NNIGHIERARCHY_HPP
