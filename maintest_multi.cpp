#include <iostream>
#include <fstream>
#include <math.h>

#include "includes.hpp"

//! \file

//! Static main program to test a multivariate hierarchy.

//! You can change the classes used for the model through the aliases below.

using HypersType = HypersFixedNNW;
using MixtureType = DirichletMixture;
template <class HypersType> using HierarchyType = HierarchyNNW<HypersType>;


int main(int argc, char *argv[]){
    std::cout << "Running maintest_multi.cpp" << std::endl;

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


    // =========================================================================
    // READ DATA AND GRID FROM FILE
    // =========================================================================
    std::string datafile = argv[1];
    Eigen::MatrixXd data = read_eigen_matrix(datafile);
    Eigen::MatrixXd grid = read_eigen_matrix("csv/grid_multi.csv");


    // =========================================================================
    // SET MODEL PARAMETERS
    // =========================================================================
    Eigen::Matrix<double,1,2> mu0;  mu0 << 5.5, 5.5;
    double lambda = 0.2;
    double nu = 5.0;
    Eigen::MatrixXd tau0 = (1/nu) * Eigen::Matrix<double, 2, 2>::Identity();
    HypersType hy(mu0, lambda, tau0, nu);

    double totalmass = 1.0;
    MixtureType mix(totalmass);


    // =========================================================================
    // LOAD ALGORITHM FACTORY
    // =========================================================================
    using Builder = std::function< std::unique_ptr<Algorithm<HierarchyType,
        HypersType, MixtureType>>(HypersType, MixtureType, Eigen::MatrixXd) >;

    Builder neal2builder = [](HypersType hy, MixtureType mix,
        Eigen::MatrixXd data){
        return std::make_unique< Neal2<HierarchyType, HypersType,
                MixtureType> >(hy, mix, data);
        };
    Builder neal8builder = [](HypersType hy, MixtureType mix,
        Eigen::MatrixXd data){
        return std::make_unique< Neal8<HierarchyType, HypersType,
                MixtureType> >(hy, mix, data);
        };

    auto &algofactory = Factory<
        Algorithm<HierarchyType, HypersType, MixtureType>, HypersType,
        MixtureType, Eigen::MatrixXd>::Instance();

    algofactory.add_builder("neal2", neal2builder);
    algofactory.add_builder("neal8", neal8builder);


    // =========================================================================
    // CREATE ALGORITHM AND SET ALGORITHM PARAMETERS
    // =========================================================================
    std::string algo = argv[2];
    auto sampler = algofactory.create_object(algo, hy, mix, data);
    (*sampler).set_rng_seed(20200229);
    (*sampler).set_maxiter(1000);
    (*sampler).set_burnin(100);


    // =========================================================================
    // CHOOSE COLLECTOR
    // =========================================================================
    BaseCollector *coll;
    std::string colltype = argv[3];
    if(colltype == "file"){
        std::string filename = "collector_multi.recordio";
        if(argc > 4){
            filename = argv[4];
        }
        else {
            std::cout << "Warning: default name " << filename <<
                " will be used for file collector" << std::endl;
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


    // =========================================================================
    // RUN SAMPLER
    // =========================================================================
    (*sampler).run(coll);


    // =========================================================================
    // DENSITY AND CLUSTER ESTIMATES
    // =========================================================================
    (*sampler).eval_density(grid, coll);
    (*sampler).write_density_to_file("csv/dens_multi.csv");
    unsigned int i_cap = (*sampler).cluster_estimate(coll);
    (*sampler).write_clustering_to_file("csv/clust_multi.csv");

    std::cout << "End of maintest_multi.cpp" << std::endl;
    return 0;
}
