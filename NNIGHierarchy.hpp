#ifndef NNIGHIERARCHY_HPP
#define NNIGHIERARCHY_HPP

#include <array>
#include <memory>
#include <random>
#include <vector>
#include <stan/math/prim/mat.hpp>

template<class Hypers> //Hypers = TupleWrapper, distro, ...
class NNIGHierarchy {
protected:
    using state_t = std::array<double,2>;

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
	std::shared_ptr<Hypers> get_hypers(){return hypers;}
    void set_state(const state_t &s){state = s;}
    void set_state(int pos, double val){state[pos] = val;}

    int get_count(){return hypers.use_count();}
	double eval_G0(double datum);
    double log_like(double datum);
	Eigen::VectorXf eval_G0(std::vector<double> datum);
    Eigen::VectorXf log_like(std::vector<double> datum);
    void draw();

    void sample_given_data(std::vector<double> data);


  std::vector<double> normal_gamma_update(std::vector<double> data,
    double mu0, double alpha0, double beta0, double lambda0);
};

#include "NNIGHierarchy_imp.hpp"

#endif // NNIGHIERARCHY_HPP
