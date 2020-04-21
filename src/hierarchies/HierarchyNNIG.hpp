#ifndef HIERARCHYNNIG_HPP
#define HIERARCHYNNIG_HPP

#include <array>
#include <memory>
#include <random>
#include <vector>
#include <stan/math/prim/mat.hpp>

// Normal likelihoood, Normal-InverseGamma hierarchy, that is:
// f ~ N(mu,sig^2)
// (mu,sig^2) ~ G
// G ~ DP(M, G0)  with G0 = N-IG


template<class Hypers>
class HierarchyNNIG {
protected:
    std::mt19937 rng;
    std::vector<Eigen::MatrixXd> state; // curr vals for F parameters: mu, sigma
    std::shared_ptr<Hypers> hypers; // current values for G0's parameters:
                                    // mu_0,Lambda0, alpha, beta

    std::vector<double> normal_gamma_update(const Eigen::VectorXd &data,
        const double mu0, const double alpha0, const double beta0,
        const double lambda0);

public:
    bool is_multivariate() const {return false;}
    
    // Contructors
    ~HierarchyNNIG() = default;

    HierarchyNNIG(std::shared_ptr<Hypers> hypers_): hypers(hypers_),
        state(2,Eigen::MatrixXd(1,1)){

    	state[0](0,0) = 0;
    	state[1](0,0) = 1;
    }

    // Getters and setters
    std::vector<Eigen::MatrixXd> get_state() const {return state;}
    std::shared_ptr<Hypers> get_hypers() const {return hypers;}
    void set_state(const std::vector<Eigen::MatrixXd> &state_){state = state_;}
    void set_state(int pos, Eigen::MatrixXd val){state[pos] = val;}

    Eigen::VectorXd eval_marg(const Eigen::MatrixXd &datum);
    Eigen::VectorXd like(const Eigen::MatrixXd &datum);

    void draw();

    void sample_given_data(const Eigen::MatrixXd &data);

};

#include "HierarchyNNIG.imp.hpp"

#endif // HIERARCHYNNIG_HPP
