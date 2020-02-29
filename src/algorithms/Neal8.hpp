#ifndef NEAL8_HPP
#define NEAL8_HPP

#include <fstream>
#include <math.h>
#include <tuple>
#include <vector>

#include <Eigen/Dense>
#include <stan/math/prim/mat.hpp>

#include "Neal2.hpp"
#include "../hyperparameters/HypersFixedNNIG.hpp"
#include "../hierarchies/HierarchyNNIG.hpp"
#include "../../output.pb.h"
#include "../mixtures/DirichletMixture.hpp"

template<template <class> class Hierarchy, class Hypers, class Mixture>
class Neal8: public Neal2<Hierarchy, Hypers, Mixture>{
private:
    // Mehtods parameters
    unsigned int n_aux = 3;

    // Data and values containers
    std::vector<Hierarchy<Hypers>> aux_unique_values;

    // Algorithm functions
    const void print_startup_message() override;

    void sample_allocations() override;

public:
    // Constructors and destructors
    ~Neal8() = default;
    Neal8(const std::vector<double> &data, int num_clusters, int n_aux,
        const Mixture &mix, const Hypers &hy):
        data(data), num_clusters(num_clusters), n_aux(n_aux), mixture(mix) {
            Hierarchy<Hypers> hierarchy(std::make_shared<Hypers> (hy));
            for(int h = 0; h < num_clusters; h++){
                unique_values.push_back(hierarchy);
            }
            for(int h = 0; h < n_aux; h++){
                aux_unique_values.push_back(hierarchy);
            }
    }

    // If no # initial clusters is given, it will be set equal to the data size
    Neal2(const std::vector<double> & data, const int num_clusters,
        const Mixture &mix, const Hypers &hy) :
        Algorithm(data, num_clusters, mix, hy) {}

    Neal8(std::vector<double> &data, int n_aux, const Mixture &mix,
        const Hypers &hy): Neal8(data, data.size(), n_aux, mix, hy) {}

    // Getters and setters
    const unsigned int get_n_aux(){return n_aux;}
    void set_n_aux(const unsigned int n_aux){n_aux = n_aux;}

}; // end of Class Neal8

#include "Neal8.imp.hpp"

#endif // NEAL8_HPP
