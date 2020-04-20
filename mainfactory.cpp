#include <iostream>
#include <fstream>

#include "includes_main.hpp"
//#include "math.h"

using HypersType = HypersFixedNNIG;
using MixtureType = DirichletMixture;
template <class HypersType> using HierarchyType = HierarchyNNIG<HypersType>;
template <typename... Args> using Builder = std::function< std::unique_ptr<
    Algorithm<HierarchyType, HypersType, MixtureType>> (Args...) >;


int main(int argc, char *argv[]){
    // Model parameters
    double mu0 = 5.0;
    double lambda = 0.1;
    double alpha0 = 2.0;
    double beta0 = 2.0;
    double totalmass = 1.0;
    unsigned int n_aux = 3;
    HypersType hy(mu0, lambda, alpha0, beta0);
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
    //using Builder = std::function< std::unique_ptr< Algorithm<HierarchyType,
      //  HypersType, MixtureType> >() >;
    Builder<Eigen::VectorXd,MixtureType,HypersType > neal2builder = [](
        Eigen::VectorXd data,MixtureType mix,HypersType hy){
            return std::make_unique< Neal2<HierarchyType,HypersType,
                MixtureType> >(hy, mix, data);
        };

    auto &algoFactory = Factory<
        Algorithm<HierarchyType, HypersType, MixtureType>,
        Eigen::VectorXd,MixtureType,HypersType >::Instance();

    algoFactory.add_builder("neal2",neal2builder);
    auto list = algoFactory.list_of_known_builders();
    for (auto &el : list){
        std::cout << el << std::endl;
    }

    // Create algorithm
    auto sampler = algoFactory.create_object("neal2", hy, mix, data);

    return 0;





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
  
    (*sampler).run(f);

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

    (*sampler).eval_density(grid, f);
    (*sampler).write_density_to_file("csv/density_ex.csv");
    unsigned int i_cap = (*sampler).cluster_estimate(f);
    std::cout << "Best clustering: at iteration " << i_cap << std::endl;
    (*sampler).write_final_clustering_to_file();
    (*sampler).write_best_clustering_to_file();

    return 0;
}
