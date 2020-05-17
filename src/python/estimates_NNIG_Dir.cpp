#include <iostream>
#include <fstream>
#include <string>

#include "../../includes.hpp"

namespace NNIGDir {
    using HypersType = HypersFixedNNIG;
    using MixtureType = DirichletMixture;
    template <class HypersType> using HierarchyType = HierarchyNNIG<HypersType>;
    using Builder = std::function< std::unique_ptr<Algorithm<HierarchyType,
        HypersType, MixtureType>>(HypersType,MixtureType, Eigen::VectorXd)>;
}


int estimates_NNIG_Dir(double mu0, double lambda, double alpha0, double beta0,
    double totalmass,
    const Eigen::VectorXd &grid, const std::string &algo,
    const std::string &filecoll_name = "collector.recordio",
    const std::string &densfile = "src/python/density.csv"){

    std::cout << "Running estimates_NNIG_Dir.cpp" << std::endl;
    using namespace NNIGDir;

    // Build model components
    HypersType hy(mu0, lambda, alpha0, beta0); // 5.0 0.1 2.0 2.0
    MixtureType mix(totalmass); // 1.0

    //Eigen::VectorXd data = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(
    //    v.data(), v.size());

    // Load algorithm factory
    Builder neal2builder_dataless = [](HypersType hy, MixtureType mix){
        return std::make_unique< Neal2<HierarchyType,HypersType,
                MixtureType> >(hy, mix);
        };
    
    Builder neal8builder_dataless = [](HypersType hy, MixtureType mix){
        return std::make_unique< Neal8<HierarchyType, HypersType,
                MixtureType> >(hy, mix);
        };

    auto &algoFactory = Factory<
        Algorithm<HierarchyType, HypersType, MixtureType>, HypersType,
        MixtureType>::Instance();

    algoFactory.add_builder("neal2_dataless",neal2builder_dataless);
    algoFactory.add_builder("neal8_dataless",neal8builder_dataless);

    // Create algorithm and set algorithm parameters
    auto sampler = algoFactory.create_object(algo, hy, mix);

    // Choose memory collector
    BaseCollector *coll = new FileCollector(filecoll_name);

    // Run algorithm
    (*sampler).eval_density(grid, coll);
    (*sampler).write_density_to_file(densfile);

    std::cout << "End of estimates_NNIG_Dir.cpp" << std::endl;
    return 0;
}
