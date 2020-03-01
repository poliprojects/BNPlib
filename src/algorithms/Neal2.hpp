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

public:
    // Other tools:
    void eval_density(const std::vector<double> grid) override;

    // Destructors and constructors:
    ~Neal2() = default;

    Neal2(const std::vector<double> &data, const int num_clusters,
        const Mixture &mixture, const Hypers &hy) :
        Algorithm<Hierarchy, Hypers, Mixture>::Algorithm(data, num_clusters,
            mixture, hy) {}

    // If no # initial clusters is given, it will be set equal to the data size:
    Neal2(const std::vector<double> &data, const Mixture &mixture,
        const Hypers &hy) :
        Algorithm<Hierarchy, Hypers, Mixture>::Algorithm(data, mixture, hy) {}

}; // end of Class Neal2

#include "Neal2.imp.hpp"

#endif // NEAL2_HPP
