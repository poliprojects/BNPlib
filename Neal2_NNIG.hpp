#ifndef NEAL2_NNIG_HPP
#define NEAL2_NNIG_HPP

#include <fstream>
#include <tuple>
#include <vector>
#include <math.h>
#include <Eigen/Dense> 
#include <stan/math/prim/mat.hpp>
#include <type_traits>
#include <math.h>
#include "NNIGHierarchy.hpp"
#include "SimpleMixture.hpp"
#include "HypersFixed.hpp"
#include "output.pb.h"

// Normal likelihoood, Normal-InverseGamma hierarchy, that is:
// f ~ N(mu,sig^2)
// (mu,sig^2) ~ G
// G ~ DP(M, G0)  with G0 = N-IG

// TODO transfer over all changes from Neal8 to here

template<template <class> class Hierarchy, class Hypers, class Mixture>
class Neal2{

private:
    unsigned int maxiter = 100;
    unsigned int burnin = 0;
    std::mt19937 rng;
    int num_clusters;
    Mixture mixture;
    std::pair< std::vector<double>, Eigen::VectorXd > density;
    ChainOutput chain;
    IterationOutput best_clust;
 

    std::vector<double> data;
    std::vector<unsigned int> allocations; // the c vector
    std::vector<Hierarchy<Hypers>> unique_values;



    void initalize();

    void step(){
        sample_allocations();
        sample_unique_values();
    }

    void sample_allocations();

    void sample_unique_values();

    void save_iteration(unsigned int iter);

    void print();

public:
    // Running tool
    void run(){
        std::cout << "Running Neal2" << std::endl;
        initalize();
        unsigned int iter = 0;
        while(iter < maxiter){
            step();    
            if(iter >= burnin){
              save_iteration(iter);
            }
            iter++;
        }
        std::cout << "Done" << std::endl;
    }


	unsigned int cluster_estimate();

    void eval_density(const std::vector<double> grid);

	const void write_final_clustering_to_file(
        std::string filename = "final_clust.csv");

    const void write_best_clustering_to_file(
        std::string filename = "best_clust.csv");

    const void write_chain_to_file(
        std::string filename = "chain.csv");

    const void write_density_to_file(
        std::string filename = "density.csv");

    // Constructors and destructors:
    ~Neal2() = default;
    Neal2(const std::vector<double> & data, const int num_clusters,
        const Mixture &mix,const Hypers &hy):
        data(data), num_clusters(num_clusters), mixture(mix) {
            Hierarchy<Hypers> hierarchy(std::make_shared<Hypers> (hy));
            for(int h = 0; h < num_clusters; h++) {
                unique_values.push_back(hierarchy);
            }
            
    }
    // If no # initial clusters is given, it will be set equal to the data size:
    Neal2(std::vector<double> &data, const Mixture & mix, const Hypers &hy):
        Neal2(data, data.size(), mix, hy) {}

}; // end of Class Neal2

#include "Neal2_NNIG.imp.hpp"

#endif // NEAL2NNIG_HPP
