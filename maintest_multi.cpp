#include <iostream>
#include <fstream>
#include <chrono>
#include "math.h" // TODO < > ?

#include "includes.hpp"

using HypersType = HypersFixedNNW;
using MixtureType = DirichletMixture;
template <class HypersType> using HierarchyType = HierarchyNNW<HypersType>;


int main(int argc, char *argv[]){
    std::cout << "Running maintest_multi.cpp" << std::endl;

    // 3D-vectorial data
    Eigen::MatrixXd data = read_eigen_matrix("csv/data_multi_2cl.csv");

    // Set model parameters
    Eigen::Matrix<double,1,2> mu0;  mu0 << 5.5, 5.5;
    double lambda = 0.2;
    double nu = 5.0;
    Eigen::MatrixXd tau0 = (1/nu) * Eigen::Matrix<double, 2, 2>::Identity();
    HypersType hy(mu0, lambda, tau0, nu);

    double totalmass = 1.0;
    MixtureType mix(totalmass);

    // Create algorithm and set algorithm parameters
    Neal8<HierarchyType, HypersType, MixtureType> sampler(hy, mix, data);
    sampler.set_rng_seed(20200229);
    sampler.set_maxiter(5000);
    sampler.set_burnin(500);

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
   	std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();    
    sampler.run(coll);
    end = std::chrono::system_clock::now();

    typedef std::chrono::duration<int, std::ratio<1, 100000000>> shakes;

    int elapsed_seconds = std::chrono::duration_cast<shakes>(
        end-start).count();
    std::cout << "time:" << elapsed_seconds << std::endl;


    // Take 3D grid from file 
    Eigen::MatrixXd grid = read_eigen_matrix("csv/grid_multi.csv");

    // Density and clustering
    sampler.eval_density(grid, coll);
    sampler.write_density_to_file("csv/dens_multi.csv");
    //unsigned int i_cap = sampler.cluster_estimate(coll);
    //sampler.write_clustering_to_file("csv/clust_multi.csv");

    std::cout << "End of maintest_multi.cpp" << std::endl;
    return 0;
}
