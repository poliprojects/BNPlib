#ifndef NEAL2_HPP
#define NEAL2_HPP

#include "Algorithm.hpp"

template<template <class> class Hierarchy, class Hypers, class Mixture>
class Neal2 : public Algorithm<Hierarchy, Hypers, Mixture> {

protected:
    const void print_startup_message() override;
    void initialize() override;
    void sample_allocations() override;
    void sample_unique_values() override;
    void sample_weights() override {return;}
    void update_hypers() override {return;}

public:
    // Other tools:
    Eigen::VectorXd density_marginal_component(Hierarchy<Hypers> &temp_hier,
        unsigned int n) override;


    // Destructors and constructors:
    ~Neal2() = default;
    Neal2()=default;
    Neal2(const Hypers &hypers_, const Mixture &mixture_,
        const Eigen::MatrixXd &data_, const unsigned int init = 0) :
        Algorithm<Hierarchy, Hypers, Mixture>::Algorithm(hypers_, mixture_,
            data_, init) {}

    Neal2(const Hypers &hypers_, const Mixture &mixture_,
        const unsigned int init = 0) :
        Algorithm<Hierarchy, Hypers, Mixture>::Algorithm(hypers_, mixture_,
             init) {}

}; // end of Class Neal2

#include "Neal2.imp.hpp"

#endif // NEAL2_HPP
