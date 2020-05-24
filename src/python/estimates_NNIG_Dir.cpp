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


int estimates_NNIG_Dir(double mu0, double lambda_, double alpha0, double beta0,
    double totalmass,
    const std::string &gridfile, const std::string &algo,
    const std::string &collfile = "collector.recordio",
    const std::string &densfile = "src/python/density.csv",
    const std::string &clustfile = "src/python/clust.csv",
    const std::string &only = ""){

    std::cout << "Running estimates_NNIG_Dir.cpp" << std::endl;
    using namespace NNIGDir;

    // Build model components
    HypersType hy(mu0, lambda_, alpha0, beta0); // 5.0 0.1 2.0 2.0
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
