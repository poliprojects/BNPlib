#ifndef HIERARCHYDUMMY_HPP
#define HIERARCHYDUMMY_HPP

#include <array>
#include <memory>
#include <random>
#include <vector>
#include <stan/math/prim/mat.hpp>

// Dummmy multivariate hierarchy for testing purposes

template<class Hypers>
class HierarchyDummy {
protected:
    std::mt19937 rng;
    std::vector<Eigen::MatrixXd> state; // mu (vector), sigma (matrix)
    std::shared_ptr<Hypers> hypers; // current values for G0's parameters:
                                    // mu0 (vector) ,lambda0 (matrix)

public:
    // Contructors
    ~HierarchyDummy() = default;
    HierarchyDummy(std::shared_ptr<Hypers> hypers): hypers(hypers){
        state.push_back(Eigen::VectorXd(2.9, 2.9, 2.9));
        state.push_back(Eigen::MatrixXd::Identity(3,3));
    }

    // Getters and setters
    std::vector<Eigen::MatrixXd>get_state(){return state;}
    std::shared_ptr<Hypers> get_hypers(){return hypers;}
    void set_state(const std::vector<Eigen::MatrixXd> &s){state = s;}
    int get_count(){return hypers.use_count();} // TODO what?

    Eigen::VectorXd eval_marg(const Eigen::MatrixXd &datum);
    Eigen::VectorXd like(const Eigen::MatrixXd &datum);

    void draw();

    std::vector<Eigen::MatrixXd> dummy_update(const Eigen::MatrixXd &data,
        const Eigen::VectorXd &mu0, const Eigen::MatrixXd &lambda0);

    void sample_given_data(const Eigen::MatrixXd &data);
};

#include "HierarchyDummy.imp.hpp"

#endif // HIERARCHYDUMMY_HPP
