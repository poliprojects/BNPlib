#include <iostream>
#include <fstream>
#include <chrono>

#include "includes.hpp"
#include "lib/math/lib/eigen_3.3.3/Eigen/Dense"

using HypersType = HypersFixedNNW;
using MixtureType = DirichletMixture;
template <class HypersType> using HierarchyType = HierarchyNNW<HypersType>;

int main(int argc, char *argv[]){
    std::cout << "Running maintest_nnws.cpp" << std::endl;

    typedef std::chrono::duration<int, std::ratio<1, 100000000>> shakes;

    // 3D-vectorial data
    Eigen::MatrixXd data = read_eigen_matrix("csv/data_multi_2cl.csv");

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
    sampler.set_rng_seed(20200229);
    sampler.set_maxiter(1000);
    sampler.set_burnin(100);

    BaseCollector *coll = new MemoryCollector();

    // Run sampler
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    sampler.run(coll);
    end = std::chrono::system_clock::now();
    int time = std::chrono::duration_cast<shakes>(end-start).count();
    std::cout << "Algo time: " << time << std::endl;

    // Take 3D grid from file
    Eigen::MatrixXd grid = read_eigen_matrix("csv/grid_multi.csv");

    // Density
    start = std::chrono::system_clock::now();
    sampler.eval_density(grid, coll);
    end = std::chrono::system_clock::now();
    sampler.write_density_to_file("csv/dens_multi.csv");
    time = std::chrono::duration_cast<shakes>(end-start).count();
    std::cout << "Dens time: " << time << std::endl;

    std::cout << "End of maintest_nnws.cpp" << std::endl;
    return 0;
}
