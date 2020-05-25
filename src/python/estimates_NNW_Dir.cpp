#include <iostream>
#include <fstream>
#include <string>

#include "../../includes.hpp"

namespace NNWDir {
    using HypersType = HypersFixedNNW;
    using MixtureType = DirichletMixture;
    template <class HypersType> using HierarchyType = HierarchyNNW<HypersType>;
    
    using BuilderDL = std::function< std::unique_ptr<Algorithm<HierarchyType,
        HypersType, MixtureType>>(HypersType, MixtureType)>; 
}


int estimates_NNW_Dir(const Eigen::Matrix<double, 1, Eigen::Dynamic>  mu0, double lambda, const Eigen::MatrixXd tau0,
    double nu, double totalmass,
    const Eigen::MatrixXd &grid, const std::string &algo,
    const std::string &collfile = "collector_multi.recordio",
    const std::string &densfile = "src/python/density_multi.csv",
    const std::string &clustfile = "src/python/clust_multi.csv",
    const std::string &only = ""){

    std::cout << "Running estimates_NNW_Dir.cpp" << std::endl;
    using namespace NNWDir;

    // Build model components
    HypersType hy(mu0, lambda, tau0, nu);
    MixtureType mix(totalmass); // 1.0

    // Load algorithm factory
    BuilderDL neal2builder_dataless = [](HypersType hy, MixtureType mix){
        return std::make_unique< Neal2<HierarchyType, HypersType,
                MixtureType> >(hy, mix);
        }; 
    BuilderDL neal8builder_dataless = [](HypersType hy, MixtureType mix){
        return std::make_unique< Neal8<HierarchyType, HypersType,
                MixtureType> >(hy, mix);
        };

    auto &algoFactory = Factory<
        Algorithm<HierarchyType, HypersType, MixtureType>, HypersType,
        MixtureType>::Instance();

    algoFactory.add_builder("neal2_dataless", neal2builder_dataless);
    algoFactory.add_builder("neal8_dataless", neal8builder_dataless);

    // Create algorithm
    auto sampler = algoFactory.create_object(algo, hy, mix);

    // Choose memory collector
    BaseCollector *coll = new FileCollector(collfile);

    // Compute estimates
    if(only != "clust"){
        (*sampler).eval_density(grid, coll);
        (*sampler).write_density_to_file(densfile);
    }
    if(only != "dens"){
        unsigned int i_cap = (*sampler).cluster_estimate(coll);
        (*sampler).write_clustering_to_file(clustfile);
    }

    std::cout << "End of estimates_NNIG_Dir.cpp" << std::endl;
    return 0;
}
