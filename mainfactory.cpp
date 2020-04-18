#include <iostream>
#include <fstream>

//#include <boost/random/random_number_generator.hpp>
//#include <boost/random/detail/qrng_base.hpp>
#include <chrono>

#include "includes_main.hpp"
#include "math.h"

int main(int argc, char *argv[]){

    load_factory();
    auto &algoFactory = Factory::Instance();
    auto list = algoFactory.registered();
    for (auto &el : list){
        std::cout << el << std::endl;
    }
    auto object = algoFactory.create("neal2");
    //object.run();

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
 
    HypersFixedNNIG hy(5.0, 1.0, 2.0, 2.0); // mu0, lambda, alpha0, beta0
    DirichletMixture mix(1); // total mass

    Neal8<HierarchyNNIG, HypersFixedNNIG, DirichletMixture> sampler(
        data, 3, mix, hy); // n_aux = 3
  
    // Run sampler
    BaseCollector *f;
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
            filename = "collector.bin";
        }
        else {
            std::string filename = argv[2];
            if(argc > 4){
                std::cout << "Warning: unused extra args present" << std::endl;
            }
        }
        f = new FileCollector(filename);
    }

    else if(collector == "memory"){
        if(argc > 3){
            std::cout << "Warning: unused extra args present" << std::endl;
        }
        f = new MemoryCollector();
    }

    else {
        std::cerr << "Error: collector type arg must be \"file\" or \"memory\""
            << std::endl;
        return 1;
    }
  
    sampler.run(f);

    // Density and clustering stuff
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

    sampler.eval_density(grid, f);
    sampler.write_density_to_file("csv/density_ex.csv");
    unsigned int i_cap = sampler.cluster_estimate(f);
    std::cout << "Best clustering: at iteration " << i_cap << std::endl;
    sampler.write_final_clustering_to_file();
    sampler.write_best_clustering_to_file();

    return 0;
}
