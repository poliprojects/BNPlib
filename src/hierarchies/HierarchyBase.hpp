#ifndef HIERARCHYBASE_HPP
#define HIERARCHYBASE_HPP

#include <array>
#include <memory>
#include <random>
#include <vector>
#include <stan/math/prim/fun.hpp>

#include <stan/math/prim/prob.hpp>
template<class Hypers>
class HierarchyBase {
protected:
    std::mt19937 rng;
    std::vector<Eigen::MatrixXd> state;
    std::shared_ptr<Hypers> hypers;

    virtual void check_state_validity() = 0;

public:
    virtual bool is_multivariate() const = 0;

    // Constructor
    HierarchyBase() = default;

    // Destructor
    virtual ~HierarchyBase() = default;

    // Getters and setters
    std::vector<Eigen::MatrixXd> get_state() const {return state;}
    std::shared_ptr<Hypers> get_hypers() const {return hypers;}
    virtual void set_state(const std::vector<Eigen::MatrixXd> &state_,
        bool check = true){
        state = state_;
        if(check){
            check_state_validity();
        }
    }

    void set_rng_seed(const unsigned int seed){rng.seed(seed);}

    virtual Eigen::VectorXd eval_marg(const Eigen::MatrixXd &data) = 0;
    virtual Eigen::VectorXd like(const Eigen::MatrixXd &data) = 0;

    virtual void draw() = 0;

    virtual void sample_given_data(const Eigen::MatrixXd &data) = 0;
};


#endif // HIERARCHYBASE_HPP
