#include <iostream>
#include <fstream>
#include <string>

#include "../../includes.hpp"




int run_NNIG_Dir(double mu0, double lambda, double alpha0, double beta0,
    double totalmass,
    std::string datafile, std::string algo, std::string coll_type,
    std::string filecoll_name = "collector.recordio",
    unsigned int rng = 0, unsigned int maxit = 0, unsigned int burn = 0){

    std::cout << "Running run_NNIG_Dir.cpp" << std::endl;

    // Set aliases
    using HypersType = HypersFixedNNIG;
    using MixtureType = DirichletMixture;
    template <class HypersType> using HierarchyType = HierarchyNNIG<HypersType>;
    using Builder = std::function< std::unique_ptr<Algorithm<HierarchyType,
        HypersType, MixtureType>>(HypersType,MixtureType, Eigen::VectorXd)>;

    // Build model components
    HypersType hy(mu0, lambda, alpha0, beta0); // 5.0 0.1 2.0 2.0
    MixtureType mix(totalmass); // 1.0

    // Read data from main arg
    std::ifstream file;
    file.open(datafile);
    if(!file.is_open()){
        std::cerr << "Error: " << datafile << " file does not exist" <<
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

    if(rng != 0)   (*sampler).set_rng_seed(rng);
    if(maxit != 0) (*sampler).set_maxiter(maxit);
    if(burn != 0)  (*sampler).set_burnin(burn);

    // Choose memory collector
    BaseCollector *coll;
    if(coll_type == "file"){
        coll = new FileCollector(filecoll_name);
    }
    else if(coll_type == "memory"){
        coll = new MemoryCollector();
    }
    else {
        std::cerr << "Error: collector type arg must be \"file\" or \"memory\""
            << std::endl;
        return 1;
    }

    // Run algorithm
    (*sampler).run(coll);

    std::cout << "End of run_NNIG_Dir.cpp" << std::endl;
    return 0;
}
