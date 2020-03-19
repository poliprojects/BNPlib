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
    void eval_density(const Eigen::MatrixXd grid) override;

    // Destructors and constructors:
    ~Neal8() = default;

    Neal8(const Eigen::MatrixXd &data, const int num_clusters,
        const int n_aux, const Mixture &mixture, const Hypers &hy) :
        Neal2<Hierarchy, Hypers, Mixture>::Neal2(data, num_clusters,
        mixture, hy), n_aux(n_aux) {

        for(unsigned int i = 0; i < n_aux; i++){
            
            aux_unique_values.push_back(this->unique_values[0]);
        }
    }

    // If no # initial clusters is given, it will be set equal to the data size
    Neal8(const Eigen::MatrixXd &data, const int n_aux,
        const Mixture &mixture, const Hypers &hy) :
        Neal8(data, data.rows(), n_aux, mixture, hy) {}

    // Getters and setters
    const unsigned int get_n_aux(){return n_aux;}
    void set_n_aux(const unsigned int n_aux){n_aux = n_aux;}

}; // end of Class Neal8

#include "Neal8.imp.hpp"

#endif // NEAL8_HPP
