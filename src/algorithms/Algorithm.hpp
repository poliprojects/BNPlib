#ifndef ALGORITHM_HPP
#define ALGORITHM_HPP

#include <fstream>
#include <math.h>
#include <random> //for std::mt19937
#include <tuple>
#include <vector>

#include <Eigen/Dense> 
#include <stan/math/prim/mat.hpp>

#include "../api/output.pb.h"
#include "../api/collectors.hpp"

template<template <class> class Hierarchy, class Hypers, class Mixture>
class Algorithm{

protected:
    // Mehtods parameters
    unsigned int maxiter = 100;
    unsigned int burnin = 50;
    int num_clusters;

    // Data and values containers
    Eigen::MatrixXd data;
    std::vector<unsigned int> allocations;
    std::vector<Hierarchy<Hypers>> unique_values;
    std::pair< Eigen::MatrixXd, Eigen::VectorXd > density; // TODO w/ Eigen
    Mixture mixture;
    ChainOutput chain;
    IterationOutput best_clust;

    // Random engine
    std::mt19937 rng;

    // Algorithm functions
    virtual const void print_startup_message() = 0;

    virtual void initialize() = 0;

    void step(){ // TODO is it virtual?
        sample_allocations();
        sample_unique_values();
        // TODO sample_weights() etc?
    }

    virtual void sample_allocations() = 0;

    virtual void sample_unique_values() = 0;

    
    IterationOutput get_state_as_proto(unsigned int iter);

    Eigen::MatrixXd proto_param_to_matrix(Param);

    const void print_ending_message();

    void save_state(BaseCollector* collector, unsigned int iter){
        collector->collect( get_state_as_proto(iter));
    }


    // Auxiliary tools
    const void print_state();    

public:
    // Running tool
    void run(BaseCollector* collector){      
        print_startup_message();
        initialize();
        unsigned int iter = 0;
        while(iter < maxiter){
            step();    
            if(iter >= burnin){
              save_state(collector,iter);
              
            }
            iter++;
        }
        print_ending_message();
    }

    //void run_and_save_cards(){ // TODO la blastiamo?
    //std::ofstream file;
    //file.open("clust_cardinalities.csv");
    //    print_startup_message();
    //    initialize();
    //    unsigned int iter = 0;
    //    while(iter < maxiter){
    //        std::cout << "Iteration # " << iter << " / " <<
    //            maxiter << std::endl; // DEBUG
    //        step();
    //        if(iter >= burnin){
    //          save_iteration(iter);
    //          file << unique_values.size() << ",";
    //        }
    //        iter++;
    //    }
    //    print_ending_message();
    //    file << std::endl;
    //    file.close();
    //}

    // Other tools
    unsigned int cluster_estimate(MemoryCollector* collector);

    virtual void eval_density(const Eigen::MatrixXd grid, MemoryCollector* collector) = 0;

    const void write_final_clustering_to_file(
        std::string filename = "clust_final.csv");

    const void write_best_clustering_to_file(
        std::string filename = "clust_best.csv");

    const void write_chain_to_file(
        std::string filename = "chain.csv");

    const void write_density_to_file(
        std::string filename = "density.csv");

    // Destructors and constructors:
    virtual ~Algorithm() = default;

    Algorithm(const Eigen::MatrixXd &data, const int num_clusters,
        const Mixture &mixture, const Hypers &hy) :
        data(data), num_clusters(num_clusters), mixture(mixture) {
            Hierarchy<Hypers> hierarchy(std::make_shared<Hypers> (hy));
            for(unsigned int i = 0; i < num_clusters; i++){
                unique_values.push_back(hierarchy);
            }
            
    }

    // If no # initial clusters is given, it will be set equal to the data size:
    Algorithm(const Eigen::MatrixXd &data, const Mixture &mixture,
        const Hypers &hy) :
        Algorithm(data, data.rows(), mixture, hy) {}

    // Getters
    const unsigned int get_maxiter(){return maxiter;}
    const unsigned int get_burnin(){return burnin;}
    const unsigned int get_num_clusters(){return num_clusters;}
    const unsigned int get_num_clusters_best(){
        return best_clust.uniquevalues_size();
    }

    // Setters
    void set_maxiter(const unsigned int maxiter){maxiter = maxiter;}
    void set_burnin(const unsigned int burnin){burnin = burnin;}
    void set_num_clusters(const unsigned int num_clusters){
        num_clusters = num_clusters;
    }

};

#include "Algorithm.imp.hpp"

#endif // ALGORITHM_HPP
