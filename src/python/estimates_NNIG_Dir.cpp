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
    std::string gridfile, std::string algo,
    std::string filecoll_name = "collector.recordio",
    std::string densfile = "src/python/density.csv"){

    std::cout << "Running estimates_NNIG_Dir.cpp" << std::endl;
    using namespace NNIGDir;

    // Build model components
    HypersType hy(mu0, lambda, alpha0, beta0); // 5.0 0.1 2.0 2.0
    MixtureType mix(totalmass); // 1.0

    // Read data from file
    std::ifstream file;
    file.open(gridfile);
    if(!file.is_open()){
        std::cerr << "Error: " << gridfile << " file does not exist" <<
            std::endl;
        return 1;
    }
    std::string str, str2;
    std::getline(file, str);
    std::istringstream iss(str);
    std::vector<double> v;
    while(std::getline(iss, str2, ',')){
        double val = ::atof(str2.c_str());
        v.push_back(val);
    }
    file.close();
    Eigen::VectorXd data = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(
        v.data(), v.size());

    // Load algorithm factory
    Builder neal2builder = [](HypersType hy, MixtureType mix,
        Eigen::VectorXd data){
        return std::make_unique< Neal2<HierarchyType,HypersType,
                MixtureType> >(hy, mix, data);
        };
    
    Builder neal8builder = [](HypersType hy, MixtureType mix,
        Eigen::VectorXd data){
        return std::make_unique< Neal8<HierarchyType,HypersType,
                MixtureType> >(hy, mix, data);
        };

    auto &algoFactory = Factory<
        Algorithm<HierarchyType, HypersType, MixtureType>, HypersType,
        MixtureType>::Instance();

    algoFactory.add_builder("neal2",neal2builder);
    algoFactory.add_builder("neal8",neal8builder);

    // Create algorithm and set algorithm parameters
    auto sampler = algoFactory.create_object(algo, hy, mix, data);

    // Choose memory collector
    BaseCollector *coll = new FileCollector(filecoll_name);

    // Run algorithm
    (*sampler).eval_density(grid, coll);
    (*sampler).write_density_to_file(densfile);

    std::cout << "End of estimates_NNIG_Dir.cpp" << std::endl;
    return 0;
}
