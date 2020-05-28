#ifndef HIERARCHYBASE_HPP
#define HIERARCHYBASE_HPP

#include <array>
#include <memory>
#include <random>
#include <vector>
#include <stan/math/prim/fun.hpp>
#include <stan/math/prim/prob.hpp>


//! Abstract base template class for a hierarchy object.

//! ...
//!
//! ...

//! \param Hypers Name of the hyperparameters class

template<class Hypers>
class HierarchyBase {
protected:
    //! Current values of the state of the Markov chain
    std::vector<Eigen::MatrixXd> state;
    //! Pointer to the hyperparameters object of the hierarchy parameters
    std::shared_ptr<Hypers> hypers;
    //! Random engine
    std::mt19937 rng;

    // AUXILIARY TOOLS
    //! Raises error if parameters are not valid w.r.t. their own domain
    virtual void check_state_validity() = 0;

public:
    //! Returns true if the hierarchy models multivariate data
    virtual bool is_multivariate() const = 0;

    // DESTRUCTOR AND CONSTRUCTORS
    virtual ~HierarchyBase() = default;
    HierarchyBase() = default;

    // EVALUATION FUNCTIONS
    //! Evaluates model likelihood in the given points
    virtual Eigen::VectorXd like(const Eigen::MatrixXd &data) = 0;
    //! Evaluates marginal distribution of data in the given points
    virtual Eigen::VectorXd eval_marg(const Eigen::MatrixXd &data) = 0;

    // SAMPLING FUNCTIONS
    //! ... (TODO check)
    virtual void draw() = 0;
    //! ...
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
