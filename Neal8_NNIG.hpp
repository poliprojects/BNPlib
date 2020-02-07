#ifndef NEAL8_NNIG_HPP
#define NEAL8_NNIG_HPP

#include <tuple>
#include <vector>
#include <Eigen/Dense>
#include <stan/math/prim/mat.hpp>
#include <type_traits>
#include "includes_universal.hpp"
#include <math.h>
#include "NNIGHierarchy.hpp"
#include "SimpleMixture.hpp"
#include "HypersFixed.hpp"
#include "output.pb.h"
#include <fstream>
// N-NIG model == gaussian kernel + N-IG base measure:
// f ~ N(mu,sig^2)
// (mu,sig^2) ~ G
// G ~ DP(M, G0)  with G0 = N-IG


// Normal likelihoood, Normal Inverse Gamma hierarchy


template<template <class> class Hierarchy, class Hypers, class Mixture>
class Neal8{
private:
    unsigned int n_aux = 3;
    unsigned int maxiter = 5000; // TODO!!!!!!!!!!!!
    unsigned int burnin = 1000;
    std::mt19937 rng;
    int numClusters;
    Mixture mixture;
	ChainOutput chain;

    std::vector<data_t> data;
    std::vector<unsigned int> allocations; // the c vector
    std::vector<Hierarchy<Hypers>> unique_values;
    std::vector<Hierarchy<Hypers>> aux_unique_values;
    std::pair< std::vector<double>, Eigen::VectorXf > density;

    void initalize();

    void step(){
        sample_allocations();
        sample_unique_values();
    }

    void sample_allocations();

    void sample_unique_values();
	
	void cluster_estimate();

    void save_iteration(unsigned int iter);

    void print();

public:
	// Running tool
    void run(){
        initalize();
        unsigned int iter = 0;
        while(iter < maxiter){
            step();
            if(iter >= burnin)
              save_iteration(iter);
            iter++;
        }
	cluster_estimate();
    }

    // Constructors and destructors:
    ~Neal8() = default;
    Neal8(const std::vector<data_t> & data, int numClusters,int n_aux,
        const Mixture & mix,const Hypers &hy):
        data(data), numClusters(numClusters), n_aux(n_aux), mixture(mix) {
            Hierarchy<Hypers> hierarchy(std::make_shared<Hypers> (hy));
            for (int h=0; h<numClusters; h++) {
                unique_values.push_back(hierarchy);
            }
            for (int h=0; h<n_aux; h++) {
                aux_unique_values.push_back(hierarchy);
            }
    }

    // If no # initial clusters is given, it will be set equal to the data size:
    Neal8(std::vector<data_t> &data, int n_aux, const Mixture & mix,
        const Hypers &hy): Neal8(data, data.size(), n_aux, mix, hy) {}

    void eval_density(const std::vector<data_t> grid);

    void write_clustering_to_file(std::string filename="output.csv");

    void write_density_to_file(std::string filename="density.csv");

}; // end of Class Neal8

#include "Neal8_NNIG_imp.hpp"

#endif // NEAL8NNIG_HPP
