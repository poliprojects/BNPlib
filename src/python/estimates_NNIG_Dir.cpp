#include <iostream>
#include <fstream>
#include <string>

#include "../../includes.hpp"

namespace NNIGDir {
    using HypersType = HypersFixedNNIG;
    using MixtureType = DirichletMixture;
    template <class HypersType> using HierarchyType = HierarchyNNIG<HypersType>;
    
    using BuilderDL = std::function< std::unique_ptr<Algorithm<HierarchyType,
        HypersType, MixtureType>>(HypersType, MixtureType)>; 
}


int estimates_NNIG_Dir(double mu0, double lambda, double alpha0, double beta0,
    double totalmass,
    const std::string &gridfile, const std::string &algo,
    const std::string &filecoll_name = "collector.recordio",
    const std::string &densfile = "src/python/density.csv",
    const std::string &clusterfile = "src/python/best_clustering.csv"){

    std::cout << "Running estimates_NNIG_Dir.cpp" << std::endl;
    using namespace NNIGDir;

    // Build model components
    HypersType hy(mu0, lambda, alpha0, beta0); // 5.0 0.1 2.0 2.0
    MixtureType mix(totalmass); // 1.0


    // Read grid from file
    std::ifstream file;
    file.open(gridfile);
    if(!file.is_open()){
        std::cerr << "Error: " << gridfile << " file does not exist" <<
            std::endl;
        return 1;
    }
    std::string str;
    std::vector<double> v;
    while(std::getline(file, str)){
        double val = ::atof(str.c_str());
        v.push_back(val);
    }
    file.close();
    Eigen::VectorXd grid = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(
        v.data(), v.size());


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

    // Create algorithm and set algorithm parameters
    auto sampler = algoFactory.create_object(algo, hy, mix);

    // Choose memory collector
    BaseCollector *coll = new FileCollector(filecoll_name);

    // Run algorithm
    (*sampler).eval_density(grid, coll);
    (*sampler).write_density_to_file(densfile);
    unsigned int i_cap = (*sampler).cluster_estimate(coll);
    (*sampler).write_clustering_to_file(clusterfile);
    
    std::cout << "End of estimates_NNIG_Dir.cpp" << std::endl;
    return 0;
}
