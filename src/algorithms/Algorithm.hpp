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

    // Data and values containers
    Eigen::MatrixXd data;
    unsigned int num_clusters;
    std::vector<unsigned int> allocations;
    std::vector< Hierarchy<Hypers> > unique_values;
    std::pair< Eigen::MatrixXd, Eigen::VectorXd > density;
    Mixture mixture;
    State best_clust;

    // Random engine
    std::mt19937 rng;

    // Flags for writing functions
    bool density_was_computed = false;
    bool clustering_was_computed = false;

    // Algorithm functions
    virtual const void print_startup_message() = 0;
    virtual void initialize() = 0;
    virtual void sample_allocations() = 0;
    virtual void sample_unique_values() = 0;
    virtual void sample_weights() = 0;
    virtual void update_hypers() = 0;
    virtual const void print_ending_message();

    void save_state(BaseCollector* collector, unsigned int iter){
        collector->collect( get_state_as_proto(iter) );
    }
    
    // Auxiliary tools
    const void print_state(); // TODO is it needed anymore?

    State get_state_as_proto(unsigned int iter);

    Eigen::MatrixXd proto_param_to_matrix(const Param& par) const;

    // Single step of algorithm
    void step(){
        sample_allocations();
        sample_unique_values();
        sample_weights();
        update_hypers();
    }

public:
    // Running tool
    void run(BaseCollector* collector){      
        print_startup_message();
        initialize();
        unsigned int iter = 0;
        while(iter < maxiter){
            step();
            if(iter >= burnin){
              save_state(collector, iter);
            }
            iter++;
        }
        collector->finish();
        print_ending_message();
    }

    // Other tools
    virtual unsigned int cluster_estimate(BaseCollector* collector);

    virtual void eval_density(const Eigen::MatrixXd &grid,
        BaseCollector* collector);

    virtual Eigen::VectorXd density_marginal_component(
        Hierarchy<Hypers> &temp_hier, unsigned int n) = 0;

    void write_clustering_to_file(
        std::string filename = "csv/clust_best.csv") const;

    void write_density_to_file(std::string filename = "csv/density.csv") const;

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

    Algorithm(const Hypers &hypers_, const Mixture &mixture_,
         const unsigned int num_clusters_ = 0) :
        mixture(mixture_), num_clusters(num_clusters_) {
        Hierarchy<Hypers> hierarchy( std::make_shared<Hypers>(hypers_) );
        unique_values.push_back(hierarchy);            
    }

    // Getters
    const unsigned int get_maxiter(){return maxiter;}
    const unsigned int get_burnin(){return burnin;}
    const unsigned int get_num_clusters(){return num_clusters;}
    const unsigned int get_num_clusters_best(){ // TODO?
        return best_clust.uniquevalues_size();
    }

    // Setters
    void set_maxiter(const unsigned int maxiter_){maxiter = maxiter_;}
    void set_burnin(const unsigned int burnin_){burnin = burnin_;}
    void set_num_clusters(const unsigned int num_clusters_){
        num_clusters = num_clusters_;
    }
    void set_rng_seed(const unsigned int seed){rng.seed(seed);}
};

#include "Algorithm.imp.hpp"

#endif // ALGORITHM_HPP
