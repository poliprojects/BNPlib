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
#include "../api/FileCollector.hpp"
#include "../api/MemoryCollector.hpp"

template<template <class> class Hierarchy, class Hypers, class Mixture>
class Algorithm{

protected:
    // Mehtods parameters

    unsigned int maxiter = 10000;
    unsigned int burnin  =  1000;

    unsigned int num_clusters;

    // Data and values containers
    Eigen::MatrixXd data;
    std::vector<unsigned int> allocations;
    std::vector<Hierarchy<Hypers>> unique_values;
    std::pair< Eigen::MatrixXd, Eigen::VectorXd > density;
    Mixture mixture;
    IterationOutput best_clust;

    // Random engine
    std::mt19937 rng;

    // Algorithm functions
    virtual const void print_startup_message() = 0;

    virtual void initialize() = 0;

    void step(){
        sample_allocations();
        sample_unique_values();
    }

    virtual void sample_allocations() = 0;

    virtual void sample_unique_values() = 0;

    
    IterationOutput get_state_as_proto(unsigned int iter);

    Eigen::MatrixXd proto_param_to_matrix(const Param& par) const;

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
            std::cout << iter << std::endl; // TODO DEBUG
            step();    
            if(iter >= burnin){
              save_state(collector, iter);    
            }
            iter++;
        }
        print_ending_message();
    }

    //void run_and_save_cards(){ // TODO la blastiamo?
    //std::ofstream file;
    //file.open("csv/clust_cardinalities.csv");
    //    print_startup_message();
    //    initialize();
    //    unsigned int iter = 0;
    //    while(iter < maxiter){
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
    unsigned int cluster_estimate(BaseCollector* collector);

    virtual void eval_density(const Eigen::MatrixXd &grid,
    	BaseCollector* collector) = 0;

    void write_final_clustering_to_file(
        std::string filename = "csv/clust_final.csv") const;

    void write_best_clustering_to_file(
        std::string filename = "csv/clust_best.csv") const;

    void write_density_to_file(
        std::string filename = "csv/density.csv") const;

    // Destructors and constructors:
    virtual ~Algorithm() = default;

    Algorithm(const Hypers &hypers_, const Mixture &mixture_,
        const Eigen::MatrixXd &data_, const unsigned int num_clusters_ = 0) :
        mixture(mixture_), data(data_), num_clusters(num_clusters_) {
        Hierarchy<Hypers> hierarchy( std::make_shared<Hypers>(hypers_) );
            if(hierarchy.is_multivariate() == false && data.cols() > 1){
            std::cout << "Warning: multivariate data supplied to " <<
               	"univariate hierarchy. The algorithm will run " <<
               	"correctly, but all data rows other than the first" <<
               	"one will be ignored" << std::endl;
            }
            if(num_clusters == 0){
                std::cout << "Warning: starting number of clusters will be " <<
                "set equal to the data size" << std::endl;
                num_clusters = data.size();
            }
            for(unsigned int i = 0; i < num_clusters; i++){
                unique_values.push_back(hierarchy);
            }
            
    }

    // Getters
    const unsigned int get_maxiter(){return maxiter;}
    const unsigned int get_burnin(){return burnin;}
    const unsigned int get_num_clusters(){return num_clusters;}
    const unsigned int get_num_clusters_best(){
        return best_clust.uniquevalues_size();
    }

    // Setters
    void set_maxiter(const unsigned int maxiter_){maxiter = maxiter_;}
    void set_burnin(const unsigned int burnin_){burnin = burnin_;}
    void set_num_clusters(const unsigned int num_clusters_){
        num_clusters = num_clusters_;
    }

};

#include "Algorithm.imp.hpp"

#endif // ALGORITHM_HPP
