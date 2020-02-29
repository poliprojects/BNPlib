#ifndef NEAL2_HPP
#define NEAL2_HPP

#include <fstream>
#include <tuple>
#include <vector>
#include <math.h>
#include <Eigen/Dense> 
#include <stan/math/prim/mat.hpp>
#include <math.h>
#include "Algorithm.hpp"
#include "../hierarchies/HierarchyNNIG.hpp"
#include "../mixtures/DirichletMixture.hpp"
#include "../hyperparameters/HypersFixedNNIG.hpp"
#include "../../output.pb.h"

// TODO transfer over all changes from Neal8 to here

template<template <class> class Hierarchy, class Hypers, class Mixture>
class Neal2 : public Algorithm<Hierarchy,Hypers,Mixture> {

private:
    void print_startup_message() override;

    void initialize() override;

    void sample_allocations() override;

    void sample_unique_values() override;

public:
    // Constructors and destructors:
    Neal2(const std::vector<double> & data, const int num_clusters,
        const Mixture &mix, const Hypers &hy) :
        Algorithm(data, num_clusters, mix, hy) {}

    // If no # initial clusters is given, it will be set equal to the data size:
    Neal2(std::vector<double> &data, const Mixture & mix, const Hypers &hy) :
        Algorithm(data, mix, hy) {}

}; // end of Class Neal2

#include "Neal2.imp.hpp"

#endif // NEAL2_HPP
