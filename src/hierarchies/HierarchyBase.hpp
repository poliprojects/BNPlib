#ifndef HIERARCHYBASE_HPP
#define HIERARCHYBASE_HPP

#include <array>
#include <memory>
#include <random>
#include <vector>
#include <stan/math/prim/fun.hpp>
#include <stan/math/prim/prob.hpp>


//! Abstract base template class for a hierarchy object.

//! This template class represents a hierarchy object in a generic iterative BNP
//! algorithm, that is, a single set of unique values with their own prior dis-
//! tribution attached to it. These values are part of the Markov chain's state
//! chain (which includes multiple hierarchies) and are simply referred to as
//! the state of the hierarchy. This object also corresponds to a single cluster
//! in the algorithm, in the sense that its state is the set of parameters for
//! the distribution of the data points that belong to it. Since the prior dis-
//! tribution for the state is often the same across multiple different hierar-
//! chies, the hyperparameters object is accessed via a shared pointer. Lastly,
//! any hierarchy that inherits from this class contains multiple ways of upda-
//! ting the state, either via prior or posterior distributions, and of evalua-
//! ting the distribution of the data, either its likelihood (whose parameters
//! are the state) or its marginal distribution.

//! \param Hypers Name of the hyperparameters class

template<class Hypers>
class HierarchyBase {
protected:
    //! Current unique values state of this cluster
    std::vector<Eigen::MatrixXd> state;
    //! Pointer to the hyperparameters object of the state
    std::shared_ptr<Hypers> hypers;
    //! Random engine
    std::mt19937 rng;

    // AUXILIARY TOOLS
    //! Raises error if the state values are not valid w.r.t. their own domain
    virtual void check_state_validity() = 0;

public:
    //! Returns true if the hierarchy models multivariate data
    virtual bool is_multivariate() const = 0;

    // DESTRUCTOR AND CONSTRUCTORS
    virtual ~HierarchyBase() = default;
    HierarchyBase() = default;

    // EVALUATION FUNCTIONS
    //! Evaluates the likelihood of data in the given points
    virtual Eigen::VectorXd like(const Eigen::MatrixXd &data) = 0;
    //! Evaluates the marginal distribution of data in the given points
    virtual Eigen::VectorXd eval_marg(const Eigen::MatrixXd &data) = 0;

    // SAMPLING FUNCTIONS
    //! Generates new values for state from the centering prior distribution
    virtual void draw() = 0;
    //! Generates new values for state from the centering posterior distribution
    virtual void sample_given_data(const Eigen::MatrixXd &data) = 0;

    // GETTERS AND SETTERS
    std::vector<Eigen::MatrixXd> get_state() const {return state;}
    std::shared_ptr<Hypers> get_hypers() const {return hypers;}
    //! \param state_ State value to set
    //! \param check  If true, a state validity check occurs after assignment
    virtual void set_state(const std::vector<Eigen::MatrixXd> &state_,
        bool check = true){
        state = state_;
        if(check){
            check_state_validity();
        }
    }
    void set_rng_seed(const unsigned int seed){rng.seed(seed);}
};


#endif // HIERARCHYBASE_HPP
