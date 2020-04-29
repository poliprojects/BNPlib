#include <iostream>
#include <fstream>

#include "includes_main.hpp"


using HypersType = HypersFixedNNIG;
using MixtureType = DirichletMixture;
template <class HypersType> using HierarchyType = HierarchyNNIG<HypersType>;

// Alias for factory

using Builder = std::function< std::unique_ptr<Algorithm<HierarchyType,
    HypersType, MixtureType>>(HypersType,MixtureType)>; 


int main(int argc, char *argv[]){
    std::cout << "Running main_factory_dataless.cpp" << std::endl;

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

    // Checks on main args
    std::ifstream file;
    if(argc < 2){
    	std::cerr << "Error: no id given for algo as arg" << std::endl;
        return 1;
   }
 

    // Load algorithm factory

    Builder neal2builder_dataless = [](HypersType hy ,MixtureType mix){

        return std::make_unique< Neal2<HierarchyType,HypersType,
                MixtureType> >(hy, mix);
        };
	
	Builder neal8builder_dataless = [](HypersType hy ,MixtureType mix){
	
        return std::make_unique< Neal8<HierarchyType,HypersType,
                MixtureType> >(hy, mix);
        };
    auto &algoFactory = Factory<
        Algorithm<HierarchyType, HypersType, MixtureType>, HypersType,
        MixtureType>::Instance();


    algoFactory.add_builder("neal2_dataless",neal2builder_dataless);
    algoFactory.add_builder("neal8_dataless",neal8builder_dataless);

    auto list = algoFactory.list_of_known_builders(); 
    std::cout << "List of known builders: ";
    for (auto &el : list){
        std::cout << el << " ";
    }
    std::cout << std::endl;

    // Create algorithm without data and set algorithm parameters
    auto sampler = algoFactory.create_object(argv[1], hy, mix);


    // Choose memory collector
    BaseCollector *coll;
    if(argc < 3){
        std::cerr << "Error: need file collector type " <<
            "(\"file\" or \"memory\") as arg" << std::endl;
        return 1;
    }

    std::string collector(argv[2]);
    if(collector == "file"){
        std::string filename;
        if(argc < 4){
            // Use default name
            filename = "collector.recordio";
        }
        else {
            std::string filename = argv[3]; 
            if(argc > 4){
                std::cout << "Warning: unused extra args present" << std::endl;
            }
        }
        coll = new FileCollector(filename);
    }

    else if(collector == "memory"){
        if(argc > 3){
            std::cout << "Warning: unused extra args present" << std::endl;
        }
        coll = new MemoryCollector();
    }

    else {
        std::cerr << "Error: collector type arg must be \"file\" or \"memory\""
            << std::endl;
        return 1;
    }



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

//////////////////////
//manca ricostruzione chains da proto da dare in pasto a eval_density TODO


    (*sampler).eval_density(grid, coll);
    (*sampler).write_density_to_file("csv/density_fact.csv");
    unsigned int i_cap = (*sampler).cluster_estimate(coll);
    (*sampler).write_clustering_to_file("csv/clust_fact.csv");

    std::cout << "End of main_factory_dataless.cpp" << std::endl;
    return 0;
}
