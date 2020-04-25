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

    Eigen::Matrix<double,1,3> mu0;  mu0 << 1.0, 1.0, 1.0;
    double lambda = 2.0;
    Eigen::MatrixXd tau0 = Eigen::Matrix<double, 3, 3>::Identity();
    double nu = 5.0;

    double totalmass = 1.0;

    HypersType hy(mu0, lambda, tau0, nu);
    MixtureType mix(totalmass); // total mass
    Neal8<HierarchyType, HypersType, MixtureType> sampler(hy, mix, data,
        data.rows());

    BaseCollector *f;
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
        f = new FileCollector(filename);
    }

    else if(collector == "memory"){
        if(argc > 2){
            std::cout << "Warning: unused extra args present" << std::endl;
        }
        f = new MemoryCollector();
    }

    else {
        std::cerr << "Error: first arg must be \"file\" or \"memory\"" <<
            std::endl;
        return 1;
    }
  
    // TODO DEBUG: reduce iterations
    sampler.set_maxiter(1000);
    sampler.set_burnin(100);

    sampler.run(f);

    sampler.eval_density(data, f); // TODO ?
    sampler.write_density_to_file();
    unsigned int i_cap = sampler.cluster_estimate(f);
    std::cout << "Best clustering: at iteration " << i_cap << std::endl;
    sampler.write_final_clustering_to_file();
    sampler.write_best_clustering_to_file();
    return 0;
}
