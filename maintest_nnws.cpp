#include <iostream>
#include <fstream>
#include <chrono>

#include "includes.hpp"

//! \file

//! Static main program to test different implementations of the NNW hierarchy.

//! This file uses the chrono library to test the efficiency of different func-
//! tions. You can change the classes used for the model through the aliases
//! below.
using HypersType = HypersFixedNNW;
using MixtureType = DirichletMixture;
template <class HypersType> using HierarchyType = HierarchyNNW<HypersType>;


int main(int argc, char *argv[]){
    std::cout << "Running maintest_nnws.cpp" << std::endl;

    // =========================================================================
    // READ DATA AND GRID FROM FILE
    // =========================================================================
    Eigen::MatrixXd data = read_eigen_matrix("csv/data_5d.csv");
    Eigen::MatrixXd grid = read_eigen_matrix("csv/grid_5d.csv");


    // =========================================================================
    // SET MODEL PARAMETERS
    // =========================================================================
    Eigen::Matrix<double,1,5> mu0; for(int i = 0; i < 5; i++) mu0[i] = 5.5;
    double lambda = 0.2;
    double nu = 5.0;
    Eigen::MatrixXd tau0 = (1/nu) * Eigen::Matrix<double, 5, 5>::Identity();
    HypersType hy(mu0, lambda, tau0, nu);

    double totalmass = 1.0;
    MixtureType mix(totalmass);    


    // =========================================================================
    // CREATE ALGORITHM AND SET ALGORITHM PARAMETERS
    // =========================================================================
    Neal2<HierarchyType, HypersType, MixtureType> sampler(hy, mix, data);
    sampler.set_rng_seed(20200229);
    sampler.set_maxiter(5000);
    sampler.set_burnin(500);


    // =========================================================================
    // INITIALIZE OBJECTS
    // =========================================================================
    BaseCollector *coll = new MemoryCollector();
    std::chrono::time_point<std::chrono::system_clock> start, end;
    using shakes = std::chrono::duration<int, std::ratio<1, 100000000>>;


    // =========================================================================
    // RUN SAMPLER AND PRINT TIME DURATION
    // =========================================================================
    start = std::chrono::system_clock::now();
    sampler.run(coll);
    end = std::chrono::system_clock::now();
    long unsigned int time = std::chrono::duration_cast<shakes>(
        end-start).count();
    std::cout << "Algo time: " << time << std::endl;


    // =========================================================================
    // RUN DENSITY ESTIMATE AND PRINT TIME DURATION
    // =========================================================================
    start = std::chrono::system_clock::now();
    sampler.eval_density(grid, coll);
    end = std::chrono::system_clock::now();
    sampler.write_density_to_file("csv/dens_5d.csv");
    time = std::chrono::duration_cast<shakes>(end-start).count();
    std::cout << "Dens time: " << time << std::endl;

    std::cout << "End of maintest_nnws.cpp" << std::endl;
    return 0;
}
