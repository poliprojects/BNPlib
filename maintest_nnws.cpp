#include <iostream>
#include <fstream>
#include <chrono>

#include "includes.hpp"
#include "lib/math/lib/eigen_3.3.3/Eigen/Dense"
//#include "math.h"

using HypersType = HypersFixedNNW;
using MixtureType = DirichletMixture;
template <class HypersType> using HierarchyType = HierarchyNNW<HypersType>;
template <class HypersType> using HierarchyType2 = HierarchyNNW2<HypersType>;

int main(){
    // 3D-vectorial data
    //Eigen::MatrixXd data = read_eigen_matrix("csv/data_multi_2cl.ssv");

    // Set model parameters
    Eigen::Matrix<double,1,2> mu0;  mu0 << 5.5, 5.5;
    double lambda = 0.2;
    double nu = 5.0;
    Eigen::MatrixXd tau0 = (1/nu) * Eigen::Matrix<double, 2, 2>::Identity();
    HypersType hy(mu0, lambda, tau0, nu);
    HypersType hy2(mu0, lambda, tau0, nu);
    std::shared_ptr<HypersType> hyptr = std::make_shared<HypersType>(hy);
    std::shared_ptr<HypersType> hyptr2 = std::make_shared<HypersType>(hy2);
    HierarchyType<HypersType> hier(hyptr);
    HierarchyType2<HypersType> hier2(hyptr2);

    Eigen::MatrixXd data = read_eigen_matrix("csv/data_multi_2cl.ssv");
    Eigen::VectorXd probas = hier.like(data);
    Eigen::VectorXd probas2 = hier2.like(data);
    std::cout << probas << std::endl;
    std::cout << probas2 << std::endl;
    std::cout << std::endl;
    std::cout << probas.sum() << std::endl;
    std::cout << probas2.sum() << std::endl;
    std::cout << std::endl;
    std::cout << probas / probas.sum() << std::endl;
    std::cout << probas2 / probas2.sum() << std::endl;
    return 0;
}

//int main(int argc, char *argv[]){
//    std::cout << "Running maintest_nnws.cpp" << std::endl;
//
//    // 3D-vectorial data
//    Eigen::MatrixXd data = read_eigen_matrix("csv/data_multi_2cl.ssv");
//
//    // Set model parameters
//    Eigen::Matrix<double,1,2> mu0;  mu0 << 5.5, 5.5;
//    double lambda = 0.2;
//    double nu = 5.0;
//    Eigen::MatrixXd tau0 = (1/nu) * Eigen::Matrix<double, 2, 2>::Identity();
//    HypersType hy(mu0, lambda, tau0, nu);
//    HypersType hy2(mu0, lambda, tau0, nu);
//
//    double totalmass = 1.0;
//    MixtureType mix(totalmass);
//
//    // Create algorithm and set algorithm parameters
//    Neal2<HierarchyType, HypersType, MixtureType> sampler(hy, mix, data);
//    Neal2<HierarchyType2, HypersType, MixtureType> sampler2(hy2, mix, data);
//    sampler.set_rng_seed(20200229);
//    sampler.set_maxiter(1000);
//    sampler.set_burnin(100);
//    sampler2.set_rng_seed(20200229);
//    sampler2.set_maxiter(1000);
//    sampler2.set_burnin(100);
//
//    BaseCollector *coll = new MemoryCollector();
//    BaseCollector *coll2 = new MemoryCollector();
//
//    // Run sampler 
//    sampler.run(coll);
//    std::cout << std::endl;
//    sampler2.run(coll2);
//
//    // Take 3D grid from file 
//    Eigen::MatrixXd grid = read_eigen_matrix("csv/grid_multi.ssv");
//
//    // Density and clustering
//    sampler.eval_density(grid, coll);
//    sampler.write_density_to_file("csv/dens_multi.csv");
//    sampler2.eval_density(grid, coll2);
//    sampler2.write_density_to_file("csv/dens_multi2.csv");
//
//    std::cout << "End of maintest_nnws.cpp" << std::endl;
//    return 0;
//}
