#include <iostream>
#include <fstream>

#include "includes_main.hpp"
#include "math.h"

int main(int argc, char *argv[]){
    // Dummy Test

    // 3D-vectorial data
    Eigen::MatrixXd data;
    data << 1.3, 0.9, 8.8, 2.0, -1.3,
            2.3, 5.1, 4.4, 0.0, -3.2,
            3.0, 4,0, 3.3, 1.1, +1.5;

    Eigen::VectorXd mu0(3.0, 3.0, 3.0);
    Eigen::MatrixXd lambda0 = 2 * Eigen::Matrix<double, 3, 3>::Identity();
    double totalmass = 1.0;
    int n_aux = 3;
    HypersDummy hy(mu0, lambda0);
    DirichletMixture mix(totalmass); // total mass
    Neal8<HierarchyDummy, HypersDummy, DirichletMixture> sampler(
        data, n_aux, mix, hy);

    BaseCollector *f;
    std::string collector(argv[2]);
    if(collector=="FileCollector"){
        std::string filename(argv[3]);
        f=new FileCollector(filename);}
    if(collector=="MemoryCollector"){
        f=new MemoryCollector();}
   
    sampler.run(f);

    return 0;
}
