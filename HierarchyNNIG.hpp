#ifndef HIERARCHYNNIG_HPP
#define HIERARCHYNNIG_HPP

#include <array>
#include <memory>
#include <random>
#include <vector>
#include <stan/math/prim/mat.hpp>

template<class Hypers>
class HierarchyNNIG {
protected:
    std::mt19937 rng;
    std::vector<double> state; // current values for F's parameters: mu, sigma
    std::shared_ptr<Hypers> hypers; // current values for G0's parameters:
                                    // mu_0,Lambda0, alpha, beta

public:
    // Contructors
    ~HierarchyNNIG() = default;

    HierarchyNNIG(std::shared_ptr<Hypers> hypers): hypers(hypers), state(2,1) {}

    // Getters and setters
    std::vector<double> get_state(){return state;}
    std::shared_ptr<Hypers> get_hypers(){return hypers;}
    void set_state(const std::vector<double> &s){state = s;}
    void set_state(int pos, double val){state[pos] = val;}
    int get_count(){return hypers.use_count();}

    // Computation tools
    double eval_marg(double datum);
    double log_like(double datum);

    Eigen::VectorXd eval_marg(std::vector<double> datum);
    Eigen::VectorXd log_like(std::vector<double> datum);

    void draw();

    void sample_given_data(std::vector<double> data);

    std::vector<double> normal_gamma_update(std::vector<double> data,
        double mu0, double alpha0, double beta0, double lambda0);

};

#include "HierarchyNNIG.imp.hpp"

#endif // HIERARCHYNNIG_HPP
