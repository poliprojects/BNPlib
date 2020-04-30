#ifndef NEAL2_IMP_HPP
#define NEAL2_IMP_HPP

#include "Neal2.hpp"

template<template <class> class Hierarchy, class Hypers, class Mixture>
const void Neal2<Hierarchy, Hypers, Mixture>::print_startup_message(){
    std::cout << "Running Neal2 algorithm..." << std::endl;
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
void Neal2<Hierarchy, Hypers, Mixture>::initialize(){
    this->cardinalities.reserve(this->data.rows());
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0, this->num_clusters-1);

    for(int h = 0; h < this->num_clusters; h++){
      this->allocations.push_back(h);
      this->cardinalities.push_back(1);
    }
    for(int j = this->num_clusters; j < this->data.rows(); j++){
        unsigned int clust = distribution(generator);
        this->allocations[j] = clust;
        this->cardinalities[clust] += 1;
    }
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
void Neal2<Hierarchy, Hypers, Mixture>::sample_allocations(){
    unsigned int k, n_unique, singleton;
    unsigned int n = this->data.rows();
    double M = this->mixture.get_totalmass();

    for(unsigned int i = 0; i < n; i++){
    	// for each data unit data[i]
    	Eigen::Matrix<double, 1, Eigen::Dynamic> datum = this->data.row(i);

        // Initialize cardinalities of unique values
        singleton = 0;
        n_unique = this->unique_values.size();

        if(this->cardinalities[ this->allocations[i] ] == 1){
        	// datum i is a singleton
            singleton = 1;
        }
        
        // Remove datum from cluster
        this->cardinalities[ this->allocations[i] ] -= 1;
        
        // Compute probabilities of clusters
        Eigen::VectorXd probas(n_unique+(1-singleton)); 

        double tot = 0.0;
        
        for(unsigned int k = 0; k < n_unique; k++){
            probas(k) = this->mixture.prob_existing_cluster(
            	this->cardinalities[k],n) *
                this->unique_values[k].like(datum)(0);

            if(singleton == 1 && k == i){
              probas(i) = this->mixture.prob_new_cluster(n, n_unique) *
                this->unique_values[0].eval_marg(datum)(0);
            }
 
            tot += probas(k);
        }

        if(singleton == 0){
            probas(n_unique) = this->mixture.prob_new_cluster(n, n_unique) *
                this->unique_values[0].eval_marg(datum)(0);
            tot += probas(n_unique);
        }

        // Normalize
        probas = probas / tot;
        
        // Draw a NEW value for ci
        unsigned int c_new = stan::math::categorical_rng(probas, this->rng) - 1;
        
        // Assign datum to the new cluster
        if(singleton == 1){
            if(c_new == this->allocations[i]){
                // case 1 of 4: SINGLETON - SINGLETON
                Eigen::VectorXd datum_col = datum;
                this->unique_values[ this->allocations[i] ].sample_given_data(
                	datum_col);
                this->cardinalities[c_new] += 1;
            }

            else{ // case 2 of 4: SINGLETON - CLUSTER
                this->unique_values.erase(
                    this->unique_values.begin() + this->allocations[i] );
                
                unsigned int c_old = this->allocations[i];
                this->allocations[i] = c_new;
                for(auto &c : this->allocations){ // relabeling
                    if(c > c_old){
                        c -= 1;
                    }
                }
                this->cardinalities[c_new] += 1;
                this->cardinalities.erase(this->cardinalities.begin() + c_old);
            } // end of else
        } // end of if(singleton == 1)

        else{ // if singleton == 0
            if(c_new == n_unique){ // case 3 of 4: NOT SINGLETON - SINGLETON
                Hierarchy<Hypers> new_unique(
                	this->unique_values[0].get_hypers());
        		Eigen::VectorXd datum_col = datum;
                new_unique.sample_given_data(datum_col);
                this->unique_values.push_back(new_unique); 
                this->allocations[i] = n_unique;
                this->cardinalities.push_back(1);
            }
            else{ // case 4 of 4: NOT SINGLETON - CLUSTER
                this->allocations[i] = c_new;
                this->cardinalities[c_new] += 1;
            }
        } // end of else
    } // end of for(int i = 0; i < n; i++) loop
} // end of sample_allocations()


template<template <class> class Hierarchy, class Hypers, class Mixture>
void Neal2<Hierarchy, Hypers, Mixture>::sample_unique_values(){
    this->num_clusters = this->unique_values.size();
    std::vector<std::vector<unsigned int>> clust_idxs(this->num_clusters);
    unsigned int n = this->allocations.size();
    for(int i = 0; i < n; i++){ // save different cluster in each row
        clust_idxs[ this->allocations[i] ].push_back(i);
    }

    for(int j = 0; j < this->num_clusters; j++){
    Eigen::MatrixXd curr_data(this->data.rows(), this->data.cols());
        int k = 0;
        for(auto &idx : clust_idxs[j]){ // TODO do a for loop with k
            curr_data.row(k) = this->data.row(idx); 
            k += 1;
        }
        curr_data.conservativeResize(k, Eigen::NoChange);
        this->unique_values[j].sample_given_data(curr_data);
    }


}



template<template <class> class Hierarchy, class Hypers, class Mixture>
Eigen::VectorXd Neal2<Hierarchy, Hypers, Mixture>::density_marginal_component(
    Hierarchy<Hypers> &temp_hier, unsigned int n){
    // Component from G0 (exploit conjugacy using explicit expression)
    double M = this->mixture.get_totalmass();
    Eigen::VectorXd dens_addendum(this->density.first.rows());    
    dens_addendum = M * temp_hier.eval_marg(this->density.first) / (M+n); 
    return dens_addendum;
}


#endif // NEAL2_IMP_HPP
