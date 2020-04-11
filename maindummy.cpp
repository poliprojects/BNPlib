#include <iostream>
#include <fstream>


#include "includes_main.hpp"
#include "math.h"

int main(int argc, char *argv[]){
    // 3D-vectorial data
    Eigen::MatrixXd data(3,5);
    fill_eigen_matrix_from_file(data,"csv/data_vec.csv");

    Eigen::VectorXd mu0(3);
    mu0 << 3.0, 3.0, 3.0;
    Eigen::MatrixXd lambda0 = 2 * Eigen::Matrix<double, 3, 3>::Identity();
    double totalmass = 1.0;
    int n_aux = 3;
    HypersDummy hy(mu0, lambda0);
    DirichletMixture mix(totalmass); // total mass
    Neal8<HierarchyDummy, HypersDummy, DirichletMixture> sampler(
        data, n_aux, mix, hy);

    BaseCollector *f;
    if(argc < 2){
        std::cerr << "Error: need at least one arg (\"file\" or \"memory\")" <<
            std::endl;
        return 1;
    }

    std::string collector(argv[1]);
    if(collector == "file"){
        std::string filename;
        if(argc < 3){
            // Use default name
            filename = "collector.bin";
        }
        else {
            std::string filename = argv[2];
            if(argc > 3){
                std::cout << "Warning: unused extra args present" << std::endl;
            }
        }
        f = new FileCollector(filename);
    }

    else if(collector == "memory"){
        if(argc > 2){
            std::cout << "Warning: unused extra args present" << std::endl;
        }
        f = new MemoryCollector();
    }

    else {
        std::cerr << "Error: first arg must be \"file\" or \"memory\"" <<
            std::endl;
        return 1;
    }
  
    sampler.run(f);

    return 0;
}
