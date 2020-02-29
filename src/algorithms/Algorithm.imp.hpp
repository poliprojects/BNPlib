#ifndef ALGORITHM_IMP_HPP
#define ALGORITHM_IMP_HPP

#include "Algorithm.hpp"


template<template <class> class Hierarchy, class Hypers, class Mixture>
void Algorithm<Hierarchy,Hypers,Mixture>::save_iteration(unsigned int iter){
 

    //std::cout << "Iteration # " << iter << " / " << maxiter-1 <<
    //    std::endl; // DEBUG
    IterationOutput iter_out;

    *iter_out.mutable_allocations() = {allocations.begin(), allocations.end()};

    for(int i = 0; i < unique_values.size(); i++){
        UniqueValues temp;
        for(auto &par : unique_values[i].get_state()){
            temp.add_params(par);
        }
        iter_out.add_phi();
        *iter_out.mutable_phi(i) = temp;
    }

    chain.add_state();
    *chain.mutable_state(iter-burnin) = iter_out;

    //print_state();
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
void Algorithm<Hierarchy,Hypers,Mixture>::print_state(){
    for (int h = 0; h < num_clusters; h++) {
        std::cout << "Parameters: ";

        for (auto c : unique_values[h].get_state()){
            std::cout << c << " " << std::endl;
        }
        std::cout << std::endl;
    }
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
void Algorithm<Hierarchy,Hypers,Mixture>::print_ending_message(){
    std::cout << "Done" << std::endl;
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
unsigned int Algorithm<Hierarchy, Hypers, Mixture>::cluster_estimate(){
    // also returns the index of the estimate in the chain object

    unsigned int niter = maxiter - burnin;
    Eigen::VectorXd errors(niter);
    int n = data.size();
    Eigen::MatrixXd tot_diss(n, n);
    tot_diss = Eigen::MatrixXd::Zero(n, n);
    std::vector<Eigen::MatrixXd> all_diss;
    IterationOutput temp;
    
    for(int h = 0; h < niter; h++){
        temp = *chain.mutable_state(h);
        Eigen::MatrixXd dissim(n, n);
        dissim = Eigen::MatrixXd::Zero(n, n);
        for(int i = 0; i < n; i++){
            for(int j = 0; j < i; j++){
                if(temp.allocations(i) == temp.allocations(j)){
                    dissim(i,j) = 1;
                }
            }
        }
    all_diss.push_back(dissim);
    tot_diss = tot_diss + dissim;
    }

    tot_diss = tot_diss / niter;

    for(int h = 0; h < niter; h++){
        // Compute error in Frobenius norm
        errors(h) = (tot_diss-all_diss[h]).norm();
    }
    
    //std::cout << errors << std::endl; // DEBUG
    std::ptrdiff_t i;
    int min_err = errors.minCoeff(&i);

    //std::ofstream file;
    //file.open("dissim_matr_mean.csv");
    //file << tot_diss;
    //file.close();

    //std::ofstream file2;
    //file2.open("dissim_matr_best.csv");
    //file2 << all_diss[i];
    //file2.close();

    best_clust = chain.state(i);
    std::cout << best_clust.phi_size() <<
        " clusters were found via least square minimization" << std::endl;
    return i;
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
void Algorithm<Hierarchy, Hypers, Mixture>::eval_density(
        const std::vector<double> grid){
    density.first = grid;

    Eigen::VectorXd dens(grid.size());
    double M = mixture.get_totalmass();
    int n = data.size();
    IterationOutput state;


    for(int iter = 0; iter < chain.state_size(); iter++){
        // for each iteration of the algorithm

        state = *chain.mutable_state(iter);
        std::vector<unsigned int> card(state.phi_size(),
            0); // TODO salviamoci ste card da qualche parte
        std::vector<double> params(state.phi(0).params_size());
        for(int j = 0; j < n; j++){
            card[ state.allocations(j) ] += 1;
        }
        Hierarchy<Hypers> temp_hier(unique_values[0].get_hypers());
        for(int h = 0; h < state.phi_size(); h++){
            for(int k = 0; k < state.phi(h).params_size(); k++){
                params[k] = state.phi(h).params(k);
            }
            temp_hier.set_state(params);

            dens += card[h] * temp_hier.like(grid) / (M+n);
        }

        // Component from G0 (exploit conjugacy: we use the explicit expression)
         dens += M * temp_hier.eval_marg(grid) / (M+n); 
    }

    // DEBUG:
    // for(int i = 0; i < grid.size(); i++)
    //     std::cout << dens(i) << " ";
    // std::cout << std::endl;

    density.second = dens / chain.state_size();

    //DEBUG:
    // for(int i = 0; i < grid.size(); i++)
    //     std::cout << density.second(i) << " ";
    // std::cout << std::endl;
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
const void Neal2<Hierarchy, Hypers, Mixture>::write_final_clustering_to_file(
        std::string filename){
    // number,datum,cluster,params1,params2,...
    
    std::ofstream file;
    file.open(filename);

    for(int i = 0; i < data.size(); i++){
        auto params = unique_values[ allocations[i] ].get_state();
        file << i << "," << data[i] << "," << allocations[i];
        for(int j = 0; j < params.size(); j++){
            file << "," << params[j];
        }
        file << std::endl;
    }
    file.close();
    std::cout << "Succesfully wrote to " << filename << std::endl;
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
const void Neal2<Hierarchy, Hypers, Mixture>::write_best_clustering_to_file(
    std::string filename){
    // number,datum,cluster,params1,params2,...

    std::ofstream file;
    file.open(filename);

    for(int i = 0; i < data.size(); i++){
        unsigned int ci = best_clust.allocations(i);
        file << i << "," << data[i] << "," << ci;
        for(int j = 0; j < best_clust.phi(ci).params_size(); j++){
            file << "," << best_clust.phi(ci).params(j);
        }
        file << std::endl;
    }
    file.close();
    std::cout << "Succesfully wrote to " << filename << std::endl;
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
const void Neal2<Hierarchy, Hypers, Mixture>::write_chain_to_file(
    std::string filename){
    // number,datum,cluster,params1,params2,...

    std::ofstream file;
    file.open(filename);

    // for each iteration of the algorithm
    for(int iter = 0; iter < chain.state_size(); iter++){
        // for each data point
        for(int i = 0; i < data.size(); i++){
            auto state_iter = chain.state(iter);
            unsigned int ci = state_iter.allocations(i);
            file << iter << "," << i << "," << data[i] << "," << ci;
            for(int j = 0; j < state_iter.phi(ci).params_size(); j++){
                file << "," << state_iter.phi(ci).params(j);
            }
            file << std::endl;
        }
    }

    file.close();
    std::cout << "Succesfully wrote to " << filename << std::endl;
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
const void Neal2<Hierarchy, Hypers, Mixture>::write_density_to_file(
    std::string filename){
    std::ofstream file;
    file.open(filename);

    for(int i = 0; i < density.first.size(); i++){
        file << density.first[i] << "," << density.second(i) << std::endl;
    }
    
    file.close();
    std::cout << "Succesfully wrote to " << filename << std::endl;
}

#endif // ALGORITHM_IMP_HPP
