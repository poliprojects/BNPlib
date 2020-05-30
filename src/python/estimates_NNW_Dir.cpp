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

//! Clustering and density estimates for an NNW + Dirichlet mixture model.

//! The Markov chain from which the estimates will be produced is contained in
//! a file collector whose filename is given. A user can opt to run only the
//! clustering or the density estimate (default case is both).

//! \param mu0, lambda_, tau0 Model parameters
//! \param nu, totalmass      More model parameters
//! \param grid               Matrix of vector points to evaluate the density in
//! \param algo               Name of the algorithm used
//! \param collfile           File collector that contains the chain
//! \param densfile           File to which density estimate will be printed
//! \param clustfile          File to which cluster estimate will be printed
//! \param only               "clust" or "dens" to not run the other

int estimates_NNW_Dir(const Eigen::Matrix<double, 1, Eigen::Dynamic> &mu0,
    const double lambda_, const Eigen::MatrixXd &tau0, const double nu,
    const double totalmass,
    const Eigen::MatrixXd &grid, const std::string &algo,
    const std::string &collfile = "collector_multi.recordio",
    const std::string &densfile = "src/python/density_multi.csv",
    const std::string &clustfile = "src/python/clust_multi.csv",
    const std::string &only = ""){

    std::cout << "Running estimates_NNW_Dir.cpp" << std::endl;
    using namespace NNWDir;

    // Build model components
    HypersType hy(mu0, lambda_, tau0, nu);
    MixtureType mix(totalmass);
    Eigen::MatrixXd data;
    
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

        algoFactory.add_builder("neal2",neal2builder);
        algoFactory.add_builder("neal8",neal8builder);
    }

    // Create algorithm
    auto sampler = algoFactory.create_object(algo, hy, mix, data);

    // Create file collector
    BaseCollector *coll = new FileCollector(collfile);

    // Compute estimates
    if(only != "clust"){
        (*sampler).eval_density(grid, coll);
        (*sampler).write_density_to_file(densfile);
    }
    if(only != "dens"){
        unsigned int i_cap = (*sampler).cluster_estimate(coll);
        (*sampler).write_clustering_to_file(clustfile);
    }

    std::cout << "End of estimates_NNW_Dir.cpp" << std::endl;
    return 0;
}
