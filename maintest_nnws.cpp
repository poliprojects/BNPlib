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

int main(int argc, char *argv[]){
    std::cout << "Running maintest_nnws.cpp" << std::endl;

    typedef std::chrono::duration<int, std::ratio<1, 100000000>> shakes; // TODO

    // 3D-vectorial data
    Eigen::MatrixXd data = read_eigen_matrix("csv/data_multi_2cl.ssv");

    // Set model parameters
    Eigen::Matrix<double,1,2> mu0;  mu0 << 5.5, 5.5;
    double lambda = 0.2;
    double nu = 5.0;
    Eigen::MatrixXd tau0 = (1/nu) * Eigen::Matrix<double, 2, 2>::Identity();
    HypersType hy(mu0, lambda, tau0, nu);
    HypersType hy2(mu0, lambda, tau0, nu);

    double totalmass = 1.0;
    MixtureType mix(totalmass);

    // Create algorithm and set algorithm parameters
    Neal2<HierarchyType, HypersType, MixtureType> sampler(hy, mix, data);
    Neal2<HierarchyType2, HypersType, MixtureType> sampler2(hy2, mix, data);
    sampler.set_rng_seed(19700101);
    sampler.set_maxiter(1000);
    sampler.set_burnin(100);
    sampler2.set_rng_seed(19700101);
    sampler2.set_maxiter(1000);
    sampler2.set_burnin(100);

    BaseCollector *coll = new MemoryCollector();
    BaseCollector *coll2 = new MemoryCollector();

    // Run sampler
    std::chrono::time_point<std::chrono::system_clock> start_1, end_1,
        start_2, end_2; // TODO TIME

    start_1 = std::chrono::system_clock::now(); // TODO TIME
    sampler.run(coll);
    end_1 = std::chrono::system_clock::now(); // TODO TIME
    std::cout << std::endl; // TODO DEBUG
    start_2 = std::chrono::system_clock::now(); // TODO TIME
    sampler2.run(coll2);
    end_2 = std::chrono::system_clock::now(); // TODO TIME

    // TODO TIME:
    int time_1 = std::chrono::duration_cast<shakes>(end_1-start_1).count();
    int time_2 = std::chrono::duration_cast<shakes>(end_2-start_2).count();
    std::cout << "Algorithm run: " << std::endl;
    std::cout << time_1 << " vs " << time_2 << std::endl;

    // Take 3D grid from file 
    Eigen::MatrixXd grid = read_eigen_matrix("csv/grid_multi.ssv");

    // Density and clustering
    start_1 = std::chrono::system_clock::now(); // TODO TIME
    sampler.eval_density(grid, coll);
    end_1 = std::chrono::system_clock::now(); // TODO TIME
    sampler.write_density_to_file("csv/dens_multi.csv");
    start_1 = std::chrono::system_clock::now(); // TODO TIME
    sampler2.eval_density(grid, coll2);
    end_1 = std::chrono::system_clock::now(); // TODO TIME
    sampler2.write_density_to_file("csv/dens_multi2.csv");

    // TODO TIME:
    time_1 = std::chrono::duration_cast<shakes>(end_1-start_1).count();
    time_2 = std::chrono::duration_cast<shakes>(end_2-start_2).count();
    std::cout << "Density estimation: " << std::endl;
    std::cout << time_1 << " vs " << time_2 << std::endl;

    std::cout << "End of maintest_nnws.cpp" << std::endl;
    return 0;
}
