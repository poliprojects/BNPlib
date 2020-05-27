#ifndef NEAL2_HPP
#define NEAL2_HPP

#include "Algorithm.hpp"


//! Template class for Neal's algorithm 2 for conjugate hierarchies

//! This class implements Neal's algorithm 2 that generates a Markov chain on
//! the clustering of the provided data.
//!
//! Using this algorithm implicitly assumes that the provided hierarchy class
//! represents a hierarchical model, i.e. a model for which the posterior dis-
//! tribution for the model parameters has the same form as their prior distri-
//! bution. Conjugacy is made use of in the computation of the estimated den-
//! sity's marginal component, since the marginal distribution for the data can
//! be expressed analytically.
//!
//! The basic idea for this algorithm is randomly drawing new allocations for
//! the data points according to weights that depend on the cardinalities of the
//! current clustering and on the mixture model used; this way, sometimes new
//! clusters are created and thus new unique values for them must be generated.
//! After that, unique values for each cluster are updated via the model's pos-
//! terior distribution, which again has a closed-form expression thanks to con-
//! jugacy.

//! \param Hierarchy Name of the hierarchy template class
//! \param Hypers    Name of the hyperparameters class
//! \param Mixture   Name of the mixture class

template<template <class> class Hierarchy, class Hypers, class Mixture>
class Neal2 : public Algorithm<Hierarchy, Hypers, Mixture> {
protected:
    using Algorithm<Hierarchy, Hypers, Mixture>::data;
    using Algorithm<Hierarchy, Hypers, Mixture>::cardinalities;
    using Algorithm<Hierarchy, Hypers, Mixture>::allocations;
    using Algorithm<Hierarchy, Hypers, Mixture>::unique_values;

    // AUXILIARY TOOLS
    //! Computes a part of the density estimation for the data
    Eigen::VectorXd density_marginal_component(Hierarchy<Hypers> &temp_hier,
        unsigned int n) override;

    // ALGORITHM FUNCTIONS
    const void print_startup_message() override;
    void initialize() override;
    void sample_allocations() override;
    void sample_unique_values() override;
    //! Empty: this algorithm does not use weights
    void sample_weights() override {return;}
    //! Empty: this algorithm does not update hyperparameters
    void update_hypers() override {return;}

public:
    // DESTRUCTOR AND CONSTRUCTORS
    ~Neal2() = default;
    Neal2() = default;
    Neal2(const Hypers &hypers_, const Mixture &mixture_,
        const Eigen::MatrixXd &data_, const unsigned int init = 0) :
        Algorithm<Hierarchy, Hypers, Mixture>::Algorithm(hypers_, mixture_,
            data_, init) {}
    Neal2(const Hypers &hypers_, const Mixture &mixture_,
        const unsigned int init = 0) :
        Algorithm<Hierarchy, Hypers, Mixture>::Algorithm(hypers_, mixture_,
             init) {}
};

#include "Neal2.imp.hpp"

#endif // NEAL2_HPP
