#include <iostream>
#include <fstream>

#include "includes_main.hpp"
#include "math.h"

using HypersType = HypersFixedNNW;
using MixtureType = DirichletMixture;
template <class HypersType> using HierarchyType = HierarchyNNW<HypersType>;


int main(int argc, char *argv[]){
    std::cout << "Running mainmulti.cpp" << std::endl;

    // 3D-vectorial data
    Eigen::MatrixXd data(10,3);
    fill_eigen_matrix_from_file(data, "csv/data_multi.ssv");

    // Set model parameters
    Eigen::Matrix<double,1,3> mu0;  mu0 << 1.0, 1.0, 1.0;
    double lambda = 2.0;
    Eigen::MatrixXd tau0 = Eigen::Matrix<double, 3, 3>::Identity();
    double nu = 5.0;
    HypersType hy(mu0, lambda, tau0, nu);

    double totalmass = 1.0;
    MixtureType mix(totalmass);

    // Create algorithm and set algorithm parameters
    Neal8<HierarchyType, HypersType, MixtureType> sampler(hy, mix, data,
        data.rows());
    sampler.set_rng_seed(20200229);
    sampler.set_maxiter(1000);
    sampler.set_burnin(100);

    BaseCollector *coll;
    if(argc < 2){
        std::cerr << "Error: need at least one arg (\"file\" or \"memory\")" <<
            std::endl;
        return 1;
    }

    std::string collector(argv[1]);
    if(collector == "file"){
        std::string filename;
        if(argc < 3){
            // Use default name
            filename = "collector.recordio";
        }
        else {
            std::string filename = argv[2];
            if(argc > 3){
                std::cout << "Warning: unused extra args present" << std::endl;
            }
        }
        coll = new FileCollector(filename);
    }

    else if(collector == "memory"){
        if(argc > 2){
            std::cout << "Warning: unused extra args present" << std::endl;
        }
        coll = new MemoryCollector();
    }

    else {
        std::cerr << "Error: first arg must be \"file\" or \"memory\"" <<
            std::endl;
        return 1;
    }

    // Run sampler
    sampler.run(coll);

    // Density and clustering
    sampler.eval_density(data, coll); // TODO ?
    sampler.write_density_to_file("csv/dens_multi.csv");
    unsigned int i_cap = sampler.cluster_estimate(coll);
    sampler.write_best_clustering_to_file("csv/clust_multi.csv");
    
    return 0;
}
