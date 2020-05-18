#include <iostream>
#include <fstream>
#include "math.h"

#include "includes.hpp"

using HypersType = HypersFixedNNIG;
using MixtureType = DirichletMixture;
template <class HypersType> using HierarchyType = HierarchyNNIG<HypersType>;


int main(int argc, char *argv[]){
    std::cout << "Running maintest_uni.cpp" << std::endl;

    // Read data from file
    if(argc < 2){
        std::cerr << "Error: no filename given for data as arg" << std::endl;
        return 1;
    }
    Eigen::VectorXd data = read_eigen_matrix(argv[1]);

    // Set model parameters
    double mu0 = 5.0;
    double lambda = 0.1;
    double alpha0 = 2.0;
    double beta0 = 2.0;
    //std::cout << "Insert mu0, lambda, alpha0, beta0 values:" << std::endl;
    //std::cin >> mu0 >> lambda >> alpha0 >> beta0;
    HypersType hy(mu0, lambda, alpha0, beta0);

    double totalmass = 1.0;
    //std::cout << "Insert total mass value:" << std::endl; 
    //std::cin >> totalmass; //1.0
    MixtureType mix(totalmass);

    // Create algorithm and set algorithm parameters
    Neal2<HierarchyType, HypersType, MixtureType> sampler(hy, mix, data);
    sampler.set_rng_seed(20200229);
    sampler.set_maxiter(5000);
    sampler.set_burnin(500);
      
    // Choose collector
    BaseCollector *coll;
    if(argc < 3){
        std::cerr << "Error: no collector type (\"file\" or \"memory\") " <<
            "given as arg" << std::endl;
        return 1;
    }

    std::string collector(argv[2]);
    if(collector == "file"){
        std::string filename;
        if(argc < 4){
            // Use default name
            filename = "collector.recordio";
        }
        else {
            std::string filename = argv[2];
            if(argc > 4){
                std::cout << "Warning: unused extra args present" << std::endl;
            }
        }
        coll = new FileCollector(filename);
    }
    else if(collector == "memory"){
        if(argc > 3){
            std::cout << "Warning: unused extra args present" << std::endl;
        }
        coll = new MemoryCollector();
    }

    else {
        std::cerr << "Error: collector type must be \"file\" or \"memory\""
            << std::endl;
        return 1;
    }

    // Run sampler
    sampler.run(coll);

    // Read grid from file
    Eigen::MatrixXd grid = read_eigen_matrix("csv/grid_uni.csv");

    // Density and clustering
	sampler.eval_density(grid, coll);
    sampler.write_density_to_file("csv/dens_test_uni.csv");
    unsigned int i_cap = sampler.cluster_estimate(coll);
    sampler.write_clustering_to_file("csv/clust_test_uni.csv");

    std::cout << "End of maintest_uni.cpp" << std::endl;
    return 0;
}
