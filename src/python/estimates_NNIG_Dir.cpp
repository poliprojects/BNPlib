#include <iostream>
#include <fstream>
#include <string>

#include "../../includes.hpp"


namespace NNIGDir {
    using HypersType = HypersFixedNNIG;
    using MixtureType = DirichletMixture;
    template <class HypersType> using HierarchyType = HierarchyNNIG<HypersType>;
    
    using BuilderDL = std::function< std::unique_ptr<Algorithm<HierarchyType,
        HypersType, MixtureType>>(HypersType, MixtureType)>;
}

//! Clustering and density estimates for an NNIG + Dirichlet mixture model.

//! The Markov chain from which the estimates will be produced is contained in
//! a file collector whose filename is given. A user can opt to run only the
//! clustering or the density estimate (default case is both).

//! \param mu0, lambda_, alpha0 Model parameters
//! \param beta0, totalmass     More model parameters
//! \param grid                 Vector of points to evaluate the density in
//! \param algo                 Name of the algorithm used
//! \param collfile             File collector that contains the chain
//! \param densfile             File to which density estimate will be printed
//! \param clustfile            File to which cluster estimate will be printed
//! \param only                 "clust" or "dens" to not run the other

int estimates_NNIG_Dir(const double mu0, const double lambda_,
	const double alpha0, const double beta0, const double totalmass,
    const Eigen::VectorXd &grid, const std::string &algo,
    const std::string &collfile = "collector.recordio",
    const std::string &densfile = "src/python/density.csv",
    const std::string &clustfile = "src/python/clust.csv",
    const std::string &only = ""){

    std::cout << "Running estimates_NNIG_Dir.cpp" << std::endl;
    using namespace NNIGDir;

    // Build model components
    HypersType hy(mu0, lambda_, alpha0, beta0);
    MixtureType mix(totalmass);

    // Load algorithm factory
    BuilderDL neal2builder_dataless = [](HypersType hy, MixtureType mix){
        return std::make_unique< Neal2<HierarchyType, HypersType,
                MixtureType> >(hy, mix);
        };
    BuilderDL neal8builder_dataless = [](HypersType hy, MixtureType mix){
        return std::make_unique< Neal8<HierarchyType, HypersType,
                MixtureType> >(hy, mix);
        };

    auto &algoFactory = Factory<
        Algorithm<HierarchyType, HypersType, MixtureType>, HypersType,
        MixtureType>::Instance();

    algoFactory.add_builder("neal2_dataless", neal2builder_dataless);
    algoFactory.add_builder("neal8_dataless", neal8builder_dataless);

    // Create algorithm
    auto sampler = algoFactory.create_object(algo, hy, mix);

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

    std::cout << "End of estimates_NNIG_Dir.cpp" << std::endl;
    return 0;
}
