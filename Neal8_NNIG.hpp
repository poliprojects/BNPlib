#ifndef NEAL8_NNIG_HPP
#define NEAL8_NNIG_HPP

#include <fstream>
#include <math.h>
#include <tuple>
#include <vector>

#include <Eigen/Dense>
#include <stan/math/prim/mat.hpp>

#include "HypersFixed.hpp"
#include "NNIGHierarchy.hpp"
#include "output.pb.h"
#include "SimpleMixture.hpp"

// Normal likelihoood, Normal-InverseGamma hierarchy, that is:
// f ~ N(mu,sig^2)
// (mu,sig^2) ~ G
// G ~ DP(M, G0)  with G0 = N-IG

template<template <class> class Hierarchy, class Hypers, class Mixture>
class Neal8{
private:
    unsigned int n_aux = 3;
    unsigned int maxiter = 5000; // TODO!!!!!!!!!!!!
    unsigned int burnin = 1000;
    std::mt19937 rng;
    int num_clusters;
    Mixture mixture;
	ChainOutput chain;

    std::vector<double> data;
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
    Neal8(const std::vector<double> & data, int num_clusters,int n_aux,
        const Mixture & mix,const Hypers &hy):
        data(data), num_clusters(num_clusters), n_aux(n_aux), mixture(mix) {
            Hierarchy<Hypers> hierarchy(std::make_shared<Hypers> (hy));
            for (int h=0; h<num_clusters; h++) {
                unique_values.push_back(hierarchy);
            }
            for (int h=0; h<n_aux; h++) {
                aux_unique_values.push_back(hierarchy);
            }
    }

    // If no # initial clusters is given, it will be set equal to the data size:
    Neal8(std::vector<double> &data, int n_aux, const Mixture & mix,
        const Hypers &hy): Neal8(data, data.size(), n_aux, mix, hy) {}

    void eval_density(const std::vector<double> grid);

    void write_clustering_to_file(std::string filename="output.csv");

    void write_density_to_file(std::string filename="density.csv");

}; // end of Class Neal8

#include "Neal8_NNIG.imp.hpp"

#endif // NEAL8NNIG_HPP
