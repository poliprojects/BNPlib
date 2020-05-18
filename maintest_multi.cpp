#include <iostream>
#include <fstream>
#include "math.h" // TODO < > ?

#include "includes.hpp"

using HypersType = HypersFixedNNW;
using MixtureType = DirichletMixture;
template <class HypersType> using HierarchyType = HierarchyNNW<HypersType>;


int main(int argc, char *argv[]){
    std::cout << "Running maintest_multi.cpp" << std::endl;

    // Read data from file
    if(argc < 2){
        std::cerr << "Error: no filename given for data as arg" << std::endl;
        return 1;
    }
    Eigen::MatrixXd data = read_eigen_matrix(argv[1]);

    // Set model parameters
    Eigen::Matrix<double,1,2> mu0;  mu0 << 5.5, 5.5;
    double lambda = 0.2;
    double nu = 5.0;
    Eigen::MatrixXd tau0 = (1/nu) * Eigen::Matrix<double, 2, 2>::Identity();
    HypersType hy(mu0, lambda, tau0, nu);

    double totalmass = 1.0;
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
    Eigen::MatrixXd grid = read_eigen_matrix("csv/grid_multi.csv");

    // Density and clustering
    sampler.eval_density(grid, coll);
    sampler.write_density_to_file("csv/dens_test_multi.csv");
    unsigned int i_cap = sampler.cluster_estimate(coll);
    sampler.write_clustering_to_file("csv/clust_test_multi.csv");

    std::cout << "End of maintest_multi.cpp" << std::endl;
    return 0;
}
