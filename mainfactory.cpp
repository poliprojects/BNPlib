#include <iostream>
#include <fstream>

#include "includes_main.hpp"


using HypersType = HypersFixedNNIG;
using MixtureType = DirichletMixture;
template <class HypersType> using HierarchyType = HierarchyNNIG<HypersType>;

// Aliases for factory
using func0 = std::function< std::unique_ptr<Algorithm<HierarchyType,
    HypersType, MixtureType>>(HypersType,MixtureType)>; // TODO names?
using func1 = std::function< std::unique_ptr<Algorithm<HierarchyType,
    HypersType, MixtureType>>(HypersType,MixtureType, Eigen::VectorXd)>;
using Builder = boost::variant<func0, func1>;


int main(int argc, char *argv[]){
    std::cout << "Running mainfactory.cpp" << std::endl;

    // Set model parameters
    double mu0, lambda, alpha0, beta0;
    //std::cout << "Insert mu0, lambda, alpha0, beta0 values:" << std::endl;
    //std::cin >> mu0 >> lambda >> alpha0 >> beta0;
    mu0 = 5.0; lambda = 0.1; alpha0 = 2.0; beta0 = 2.0;
    HypersType hy(mu0, lambda, alpha0, beta0);

    double totalmass;
    //std::cout << "Insert total mass value:" << std::endl; 
    //std::cin >> totalmass; //1.0
    totalmass = 1.0;
    MixtureType mix(totalmass);

    // Read data from main arg
    std::ifstream file;
    if(argc < 2){
        std::cerr << "Error: no filename given for data as arg" << std::endl;
        return 1;
    }
    file.open(argv[1]);
    if(!file.is_open()){
        std::cerr << "Error: " << argv[1] << " file does not exist" <<
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
        v.data(), v.size()); // TODO: meglio con conservative resize?

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

    Builder neal2builder_density = [](HypersType hy ,MixtureType mix){
        // che poi andrebbe nel main SOLO per eval density
        return std::make_unique< Neal2<HierarchyType,HypersType,
                MixtureType> >(hy, mix);
        };
	
	Builder neal8builder_density = [](HypersType hy ,MixtureType mix){
	
        return std::make_unique< Neal8<HierarchyType,HypersType,
                MixtureType> >(hy, mix);
        };
    auto &algoFactory = Factory<
        Algorithm<HierarchyType, HypersType, MixtureType>, HypersType,
        MixtureType>::Instance();

    algoFactory.add_builder("neal2",neal2builder);
    algoFactory.add_builder("neal8",neal8builder);
    algoFactory.add_builder("neal2density",neal2builder_density);
    algoFactory.add_builder("neal8density",neal8builder_density);

    auto list = algoFactory.list_of_known_builders();
    for (auto &el : list){
        std::cout << el << std::endl;
    }

    // Create algorithm and set algorithm parameters
    auto sampler = algoFactory.create_object(argv[2], hy, mix, data);

    sampler.set_rng_seed(20200229);
    sampler.set_maxiter(1000);
    sampler.set_burnin(100);

    return 0; // TODO DEBUG

    // Choose memory collector
    BaseCollector *coll;
    if(argc < 4){
        std::cerr << "Error: need file collector type " <<
            "(\"file\" or \"memory\") as arg" << std::endl;
        return 1;
    }

    std::string collector(argv[3]);
    if(collector == "file"){
        std::string filename;
        if(argc < 5){
            // Use default name
            filename = "collector.recordio";
        }
        else {
            std::string filename = argv[4]; 
            if(argc > 5){
                std::cout << "Warning: unused extra args present" << std::endl;
            }
        }
        coll = new FileCollector(filename);
    }

    else if(collector == "memory"){
        if(argc > 4){
            std::cout << "Warning: unused extra args present" << std::endl;
        }
        coll = new MemoryCollector();
    }

    else {
        std::cerr << "Error: collector type arg must be \"file\" or \"memory\""
            << std::endl;
        return 1;
    }

    // Run algorithm
    (*sampler).run(coll);

    // Density and clustering
    double temp = 0.0;
    double step = 0.05;
    double upp_bnd = 10.0;
    std::vector<double> v_temp;
    while(temp <= upp_bnd){
        v_temp.push_back(temp);
        temp += step;
    }
    Eigen::VectorXd grid = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(
        v_temp.data(), v_temp.size()); 

    (*sampler).eval_density(grid, coll);
    (*sampler).write_density_to_file("csv/density_fact.csv");
    unsigned int i_cap = (*sampler).cluster_estimate(coll);
    (*sampler).write_best_clustering_to_file("csv/clust_fact.csv");

    return 0;
}
