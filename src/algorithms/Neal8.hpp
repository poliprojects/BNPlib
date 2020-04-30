#ifndef NEAL8_HPP
#define NEAL8_HPP

#include "Neal2.hpp"

template<template <class> class Hierarchy, class Hypers, class Mixture>
class Neal8 : public Neal2<Hierarchy, Hypers, Mixture>{
protected:
    // Mehtods parameters
    unsigned int n_aux = 3;

    // Data and values containers
    std::vector<Hierarchy<Hypers>> aux_unique_values;

    // Algorithm functions
    const void print_startup_message() override;

    void sample_allocations() override;

public:
    // Other tools:

    Eigen::VectorXd density_marginal_component(Hierarchy<Hypers> &temp_hier,
        unsigned int n) override;

    // Destructors and constructors:
    ~Neal8() = default;
    Neal8()=default;
    Neal8(const Hypers &hypers_, const Mixture &mixture_,
        const Eigen::MatrixXd &data_, const unsigned int init = 0) :
        Neal2<Hierarchy, Hypers, Mixture>::Neal2(hypers_, mixture_, data_,
        init) {
        for(unsigned int i = 0; i < n_aux; i++){
            aux_unique_values.push_back(this->unique_values[0]);
        }
    }

    Neal8(const Hypers &hypers_, const Mixture &mixture_,
        const unsigned int init = 0) :
        Neal2<Hierarchy, Hypers, Mixture>::Neal2(hypers_, mixture_, init) {}

    // Getters and setters
    const unsigned int get_n_aux(){return n_aux;}
    void set_n_aux(const unsigned int n_aux_){n_aux = n_aux_;}

}; // end of Class Neal8

#include "Neal8.imp.hpp"

#endif // NEAL8_HPP
