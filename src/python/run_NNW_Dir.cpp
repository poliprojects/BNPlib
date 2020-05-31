#include <iostream>
#include <fstream>
#include <string>

#include "../../includes.hpp"


namespace NNWDir {
    using HypersType = HypersFixedNNW;
    using MixtureType = DirichletMixture;
    template <class HypersType> using HierarchyType = HierarchyNNW<HypersType>;

    using Builder = std::function< std::unique_ptr<Algorithm<HierarchyType,
        HypersType, MixtureType>>(HypersType,MixtureType, Eigen::MatrixXd)>;
}

//! Runs algorithm for an NNW + Dirichlet mixture model.

//! The Markov chain produced by the algorithm will be saved to a collector. A
//! file collector must be used if one wants to follow up the algorithm run with
//! the estimates run using the Python interface of this library.

//! \param mu0, lambda_, tau0 Model parameters
//! \param nu, totalmass      More model parameters
//! \param datafile           csv file from which the model data is read
//! \param algo               Name of the algorithm used
//! \param collfile           Name of file collector to which save the chain
//! \param rng                RNG seed for the algorithm
//! \param maxit              Number of iterations of the algorithm
//! \param burn               Number of burn-in (discarded) iterations
//! \param n_aux              If algo="neal8", number of auxiliary blocks used

int run_NNW_Dir(const Eigen::Matrix<double, 1, Eigen::Dynamic> &mu0,
    const double lambda_, const Eigen::MatrixXd &tau0, const double nu,
    const double totalmass,
    const std::string &datafile, const std::string &algo,
    const std::string &collfile = "collector.recordio",
    const unsigned int rng = 0, const unsigned int maxit = 0,
    const unsigned int burn = 0, const unsigned int n_aux = 0){

    std::cout << "Running run_NNW_Dir.cpp" << std::endl;
    using namespace NNWDir;

    // Build model components
    HypersType hy(mu0, lambda_, tau0, nu);
    MixtureType mix(totalmass);

    // Read data from file
    Eigen::MatrixXd data = read_eigen_matrix(datafile);

    // Load algorithm factory
    auto &algoFactory = Factory<
        Algorithm<HierarchyType, HypersType, MixtureType>, HypersType,
        MixtureType,Eigen::MatrixXd>::Instance();

    if (!algoFactory.check_existence(algo)){

        Builder neal2builder = [](HypersType hy, MixtureType mix,
            Eigen::MatrixXd data){
            return std::make_unique< Neal2<HierarchyType,HypersType,
                    MixtureType> >(hy, mix, data);
            };

        Builder neal8builder = [](HypersType hy, MixtureType mix,
            Eigen::MatrixXd data){
            return std::make_unique< Neal8<HierarchyType,HypersType,
                    MixtureType> >(hy, mix, data);
            };

        algoFactory.add_builder("neal2", neal2builder);
        algoFactory.add_builder("neal8", neal8builder);
    }

    // Create algorithm and set algorithm parameters
    auto sampler = algoFactory.create_object(algo, hy, mix, data);

    if(rng != 0)   (*sampler).set_rng_seed(rng);
    if(maxit != 0) (*sampler).set_maxiter(maxit);
    if(burn != 0)  (*sampler).set_burnin(burn);
    if(algo == "neal8" && n_aux != 0) (*sampler).set_n_aux(n_aux);

    // Create file collector
    BaseCollector *coll = new FileCollector(collfile);

    // Run algorithm
    (*sampler).run(coll);

    std::cout << "End of run_NNW_Dir.cpp" << std::endl;
    return 0;
}
