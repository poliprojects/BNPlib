#ifndef HIERARCHYBASE_HPP
#define HIERARCHYBASE_HPP

#include <array>
#include <memory>
#include <random>
#include <vector>
#include <stan/math/prim/mat.hpp>

template<class Hypers>
class HierarchyBase {
protected:
    std::mt19937 rng;
    std::vector<Eigen::MatrixXd> state; // mu (vector), sigma (matrix)
    std::shared_ptr<Hypers> hypers; // current values for G0's parameters:
                                    // mu0 (vector) ,lambda0 (matrix)

    std::vector<Eigen::MatrixXd> dummy_update(const Eigen::MatrixXd &data,
        const Eigen::VectorXd &mu0, const Eigen::MatrixXd &lambda0);

public:
    virtual bool is_multivariate() const = 0; // TODO funzione o variabile?

    // Destructor
    virtual ~HierarchyBase() = default;

    // Getters and setters
    std::vector<Eigen::MatrixXd> get_state() const {return state;}
    std::shared_ptr<Hypers> get_hypers() const {return hypers;}
    void set_state(const std::vector<Eigen::MatrixXd> &state_){state = state_;}
    void set_rng_seed(const unsigned int seed){rng.seed(seed);}

    virtual Eigen::VectorXd eval_marg(const Eigen::MatrixXd &data) = 0;
    virtual Eigen::VectorXd like(const Eigen::MatrixXd &data) = 0;

    virtual void draw() = 0;

    virtual void sample_given_data(const Eigen::MatrixXd &data) = 0;

    // TODO add virtual check_state_validity() ?
};


#endif // HIERARCHYBASE_HPP
