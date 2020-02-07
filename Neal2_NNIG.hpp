#ifndef NEAL2_NNIG_HPP
#define NEAL2_NNIG_HPP

#include <tuple>
#include <vector>
#include <math.h>
#include <Eigen/Dense> 
#include <stan/math/prim/mat.hpp>
#include <type_traits>
#include "includes_universal.hpp"
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
	//ChainOutput chain;
 

    std::vector<data_t> data;
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
        initalize();
        unsigned int iter = 0;
        while(iter < maxiter){
            step();    
            if(iter >= burnin)
              save_iteration(iter);
            iter++;
        }
    }

    // Constructors and destructors:
    ~Neal2() = default;
    Neal2(const std::vector<data_t> & data, int num_clusters,
        const Mixture & mix,const Hypers &hy):
        data(data), num_clusters(num_clusters), mixture(mix) {
            Hierarchy<Hypers> hierarchy(std::make_shared<Hypers> (hy));
            for (int h=0; h<num_clusters; h++) {
                unique_values.push_back(hierarchy);
            }
            
    }
    // If no # initial clusters is given, it will be set equal to the data size:
    Neal2(std::vector<data_t> &data, const Mixture & mix,
        const Hypers &hy): Neal2(data, data.size(), mix, hy) {}

    

}; // end of Class Neal2

#include "Neal2_NNIG_imp.hpp"

#endif // NEAL2NNIG_HPP
