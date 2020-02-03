#ifndef NNIGHIERARCHY_HPP
#define NNIGHIERARCHY_HPP

#include <array>
#include <random>
#include <vector>
#include <memory>
#include <stan/math/prim/mat.hpp>

#include "includes_universal.hpp"


template<class Hypers> //Hypers = TupleWrapper, distro, ...
class NNIGHierarchy {
protected:
    using state_t = std::array<par_t,2>;

    std::mt19937 rng;
    state_t state; // current values for F's parameters: mu, sigma
    std::shared_ptr<Hypers> hypers; // current values for G0's parameters:
                                    // mu_0,Lambda0, alpha, beta

public:
    // Contructors:
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


  parvec_t normal_gamma_update(std::vector<data_t> data,
    double mu0, double alpha0, double beta0, double lambda0);
};

#include "NNIGHierarchy_imp.hpp"

#endif // NNIGHIERARCHY_HPP
