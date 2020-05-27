#include <iostream>
#include <fstream>
#include <math.h>

#include "includes.hpp"

using HypersType = HypersFixedNNIG;
using MixtureType = DirichletMixture;
template <class HypersType> using HierarchyType = HierarchyNNIG<HypersType>;

using Builder = std::function< std::unique_ptr<Algorithm<HierarchyType,
    HypersType, MixtureType>>(HypersType,MixtureType, Eigen::VectorXd)>;

//! \file

//! main program to test an univariate hierarchy.

//! Longer explaination of the same thing.

int main(int argc, char *argv[]){
    std::cout << "Running maintest_uni.cpp" << std::endl;

    // =========================================================================
    // CHECK MAIN ARGS
    // =========================================================================
    // [0]main [1]data [2]algo [3]coll [4]filecollname
    switch(argc){
        case 1:
            std::cerr << "Error: no data filename given as arg" << std::endl;
            return 1;
        case 2:
            std::cerr << "Error: no algorithm id given as arg" << std::endl;
            return 1;
        case 3:
           std::cerr << "Error: no collector type (\"file\" or \"memory\") " <<
                "given as arg" << std::endl;
            return 1;
        default:
            break;
    }

    // Read data from file
    std::string datafile = argv[1];
    Eigen::VectorXd data = read_eigen_matrix(datafile);

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

    // Load algorithm factory
    Builder neal2builder = [](HypersType hy, MixtureType mix,
        Eigen::VectorXd data){
        return std::make_unique< Neal2<HierarchyType, HypersType,
                MixtureType> >(hy, mix, data);
        };
    
    Builder neal8builder = [](HypersType hy, MixtureType mix,
        Eigen::VectorXd data){
        return std::make_unique< Neal8<HierarchyType, HypersType,
                MixtureType> >(hy, mix, data);
        };

    auto &algofactory = Factory<
        Algorithm<HierarchyType, HypersType, MixtureType>, HypersType,
        MixtureType>::Instance();

    algofactory.add_builder("neal2", neal2builder);
    algofactory.add_builder("neal8", neal8builder);

    // Create algorithm and set algorithm parameters
    std::string algo = argv[2];
    auto sampler = algofactory.create_object(algo, hy, mix, data);
    (*sampler).set_rng_seed(20200229);
    (*sampler).set_maxiter(5000);
    (*sampler).set_burnin(500);
    
    // Choose collector
    BaseCollector *coll;
    std::string colltype = argv[3];
    if(colltype == "file"){
        std::string filename = "collector.recordio";
        if(argc > 4){
            filename = argv[4];
        }
        coll = new FileCollector(filename);
    }
    else if(colltype == "memory"){
        coll = new MemoryCollector();
    }
    else {
        std::cerr << "Error: collector type must be \"file\" or \"memory\""
            << std::endl;
        return 1;
    }

    // Run sampler
    (*sampler).run(coll);

    // Read grid from file
    Eigen::MatrixXd grid = read_eigen_matrix("csv/grid_uni.csv");

    // Density and clustering
	(*sampler).eval_density(grid, coll);
    (*sampler).write_density_to_file("csv/dens_uni.csv");
    unsigned int i_cap = (*sampler).cluster_estimate(coll);
    (*sampler).write_clustering_to_file("csv/clust_uni.csv");

    std::cout << "End of maintest_uni.cpp" << std::endl;
    return 0;
}
