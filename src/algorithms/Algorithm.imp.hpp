#ifndef ALGORITHM_IMP_HPP
#define ALGORITHM_IMP_HPP

#include "Algorithm.hpp"


// Algorithm functions

template<template <class> class Hierarchy, class Hypers, class Mixture>
const void Algorithm<Hierarchy, Hypers, Mixture>::print_ending_message(){
    std::cout << "Done" << std::endl;
}


// Auxiliary tools

template<template <class> class Hierarchy, class Hypers, class Mixture>
State Algorithm<Hierarchy, Hypers, Mixture>::get_state_as_proto(
    unsigned int iter){

    State iter_out;
    *iter_out.mutable_allocations() = {allocations.begin(), allocations.end()};

    for(size_t i = 0; i < unique_values.size(); i++){
        UniqueValues uniquevalues_temp;
        for(size_t k = 0; k < unique_values[i].get_state().size(); k++){
            Eigen::MatrixXd par_temp = unique_values[i].get_state()[k];
            Param par_temp_proto;
            for(size_t j = 0; j < par_temp.cols(); j++){
                Par_Col col_temp;
                for(size_t h = 0; h < par_temp.rows(); h++){
                    col_temp.add_elems(par_temp(h,j));
                }
                par_temp_proto.add_par_cols();
                *par_temp_proto.mutable_par_cols(j) = col_temp;     
            }
            uniquevalues_temp.add_params();
            *uniquevalues_temp.mutable_params(k) = par_temp_proto;
        }
        iter_out.add_uniquevalues();
        *iter_out.mutable_uniquevalues(i) = uniquevalues_temp;
    }
    return iter_out;
}



template<template <class> class Hierarchy, class Hypers, class Mixture>
Eigen::MatrixXd Algorithm<Hierarchy, Hypers, Mixture>::proto_param_to_matrix(
    const Param &par) const {
    Eigen::MatrixXd par_matrix = Eigen::MatrixXd::Zero(
        par.par_cols(0).elems_size(), par.par_cols_size() );
    for(size_t h = 0; h < par.par_cols_size(); h++){
        for(size_t j = 0; j < par.par_cols(h).elems_size(); j++){
            par_matrix(j,h) = par.par_cols(h).elems(j);
        }
    }
    return par_matrix;
}


// Other tools

template<template <class> class Hierarchy, class Hypers, class Mixture>
void Algorithm<Hierarchy, Hypers, Mixture>::eval_density(
    const Eigen::MatrixXd &grid, BaseCollector* collector){

    density.first = grid;
    Eigen::VectorXd dens(Eigen::MatrixXd::Zero(grid.rows(),1));
    double M = mixture.get_totalmass();
    unsigned int n;
    State state;
    std::cout << dens << std::endl; // TODO DEBUG
    for(size_t iter = 0; iter < collector->get_size(); iter++){
        // for each iteration of the algorithm
        state = collector->get_next_state();
        std::vector<unsigned int> card(state.uniquevalues_size(), 0);
        std::vector<Eigen::MatrixXd> params(
            state.uniquevalues(0).params_size() );
        if (iter == 0){
            n = state.allocations_size();
            
        }
        for(size_t j = 0; j < n; j++){
            card[ state.allocations(j) ] += 1;
        }
        std::cout << "iter" << std::endl; // TODO DEBUG
        std::cout << iter << std::endl; // TODO DEBUG
        Hierarchy<Hypers> temp_hier(unique_values[0].get_hypers());
        for(size_t h = 0; h < state.uniquevalues_size(); h++){
            for(size_t k = 0; k < state.uniquevalues(h).params_size(); k++){
                params[k] = proto_param_to_matrix(
                    state.uniquevalues(h).params(k) );
            }
            temp_hier.set_state(params, false);
   
            dens += card[h]* temp_hier.like(grid) / (M+n);
			std::cout << "dens: " << dens << std::endl; // TODO
        }
            
        dens += density_marginal_component(temp_hier,n);
                       
    }

    density.second = dens / collector->get_size();
    density_was_computed = true;
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
unsigned int Algorithm<Hierarchy, Hypers, Mixture>::cluster_estimate(
    BaseCollector* collector){
    // also returns the index of the estimate in the chain object

    unsigned int niter = maxiter - burnin;
    Eigen::VectorXd errors(niter);
    unsigned int n = data.rows();
    Eigen::MatrixXd tot_diss(n, n);
    tot_diss = Eigen::MatrixXd::Zero(n, n);
    std::vector<Eigen::MatrixXd> all_diss;
    State temp;
    
    for(size_t h = 0; h < niter; h++){
        temp = collector->get_next_state();
        Eigen::MatrixXd dissim(n, n);
        dissim = Eigen::MatrixXd::Zero(n, n);
        for(size_t i = 0; i < n; i++){
            for(size_t j = 0; j < i; j++){
                if(temp.allocations(i) == temp.allocations(j)){
                    dissim(i,j) = 1;
                }
            }
        }
    all_diss.push_back(dissim);
    tot_diss = tot_diss + dissim;
    }

    tot_diss = tot_diss / niter;

    for(size_t h = 0; h < niter; h++){
        // Compute error in Frobenius norm
        errors(h) = (tot_diss-all_diss[h]).norm();
    }
    
    std::ptrdiff_t i;
    unsigned int min_err = errors.minCoeff(&i);

    best_clust = collector->get_state(i);
    std::cout << "Optimal clustering: at iteration " << i << " with " <<
        best_clust.uniquevalues_size() << " clusters" << std::endl;

    clustering_was_computed = true;

    return i;
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
void Algorithm<Hierarchy, Hypers, Mixture>::write_clustering_to_file(
    std::string filename) const {
    // number,datum,cluster,params1,params2,...

    if(!clustering_was_computed){
        std::cerr << "Error: cannot write clustering to file; " <<
            "cluster_estimate() must be called first" << std::endl;
        return;
    }

    std::ofstream file;
    file.open(filename);

    for(size_t i = 0; i < data.rows(); i++){
        unsigned int ci = best_clust.allocations(i);
        file << i << "," << data.row(i) << "," << ci;
        for(size_t j = 0; j < best_clust.uniquevalues(ci).params_size(); j++){
            file << "," << proto_param_to_matrix( 
                best_clust.uniquevalues(ci).params(j) );
        }
        file << std::endl;
    }
    file.close();
    std::cout << "Successfully wrote clustering to " << filename << std::endl;
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
void Algorithm<Hierarchy, Hypers, Mixture>::write_density_to_file(
    std::string filename) const {

    if(!density_was_computed){
        std::cerr << "Error: cannot write density to file; eval_density() " <<
            "must be called first" << std::endl;
        return;
    }

    std::ofstream file;
    file.open(filename);

    for(size_t i = 0; i < density.first.rows(); i++){
        Eigen::VectorXd point = density.first.row(i);
        for(size_t j = 0; j < point.size(); j++){
            file << point(j) << ",";
        }
        file << density.second(i) << std::endl;
    }
    
    file.close();
    std::cout << "Successfully wrote density to " << filename << std::endl;
}


#endif // ALGORITHM_IMP_HPP
