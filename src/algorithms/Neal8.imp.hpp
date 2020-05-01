#ifndef NEAL8_IMP_HPP
#define NEAL8_IMP_HPP

#include "Neal8.hpp"


// Algorithm functions
template<template <class> class Hierarchy, class Hypers, class Mixture>
const void Neal8<Hierarchy, Hypers, Mixture>::print_startup_message(){
    std::cout << "Running Neal8 algorithm (with m=" << n_aux <<
        " auxiliary blocks)..." << std::endl;
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
void Neal8<Hierarchy, Hypers, Mixture>::sample_allocations(){
    // Initialize some relevant variables
    unsigned int n_unique, singleton;
    unsigned int n = data.rows();

    for(size_t i = 0; i < n; i++){ // for each data unit data[i]
    	Eigen::Matrix<double, 1, Eigen::Dynamic> datum = data.row(i);

        singleton = 0;
        n_unique = unique_values.size();
     
        if(cardinalities[ allocations[i] ] == 1){
        	// datum i is a singleton
            aux_unique_values[0].set_state( unique_values[
                allocations[i] ].get_state() ); // move phi value in aux
            singleton = 1;
        }

        cardinalities[ allocations[i] ] -= 1;
        
        // Draw the aux from G0
        for(size_t j = singleton; j < n_aux; j++){
            aux_unique_values[j].draw();
        }

        // Draw a NEW value for ci
        Eigen::VectorXd probas(n_unique+n_aux); //k or n_unique
        
        double tot = 0.0;
        for(size_t k = 0; k < n_unique; k++){ // if datum i is a singleton, then
            // card[k] when k=allocations[i] is equal to 0 -> probas[k]=0
            probas(k) = this->mixture.prob_existing_cluster(
            	cardinalities[k], n) * unique_values[k].like(datum)(0);
            tot += probas(k);
        }

        for(size_t k = 0; k < n_aux; k++){
            probas(n_unique+k) = this->mixture.prob_new_cluster(n, n_unique) *
                aux_unique_values[k].like(datum)(0) / n_aux;
            tot += probas(n_unique+k);
        }

        // Normalize
        probas = probas / tot;

        // Draw a NEW value for ci
        unsigned int c_new = stan::math::categorical_rng(probas, this->rng) - 1;

        if(singleton == 1){
            if(c_new >= n_unique){ // case 1 of 4: SINGLETON - AUX
                unique_values[ allocations[i] ].set_state(
                    aux_unique_values[c_new-n_unique].get_state());
                cardinalities[ allocations[i] ] += 1;
            }
            else{ // case 2 of 4: SINGLETON - OLD VALUE
                unique_values.erase( unique_values.begin() + allocations[i] );
                
                unsigned int c_old = allocations[i];
                allocations[i] = c_new;
                for(auto &c : allocations){ // relabeling
                    if(c > c_old){
                        c -= 1;
                    }
                }
                cardinalities[c_new] += 1;
                cardinalities.erase(cardinalities.begin() + c_old);
            } // end of else
        } // end of if(singleton == 1)
        
        else{ // if singleton == 0
            if(c_new >= n_unique){ // case 3 of 4: NOT SINGLETON - AUX
                unique_values.push_back( aux_unique_values[c_new-n_unique] );
                cardinalities.push_back(1);
                allocations[i] = n_unique;
            }
            else{ // case 4 of 4: NOT SINGLETON - OLD VALUES
                allocations[i] = c_new;
                cardinalities[c_new] += 1;
            }
        } // end of else
    } // end of datum[i] loop
} // end of sample_allocations()


// Other tools
template<template <class> class Hierarchy, class Hypers, class Mixture>
Eigen::VectorXd Neal8<Hierarchy, Hypers, Mixture>::density_marginal_component(
    Hierarchy<Hypers> &temp_hier, unsigned int n){
    double M = this->mixture.get_totalmass();
    Eigen::VectorXd dens_addendum(this->density.first.rows());
            
    // Component from G0
    for(size_t h = 0; h < n_aux; h++){
        temp_hier.draw();
        dens_addendum += temp_hier.like(this->density.first) * (M/n_aux)/(M+n);
    }

    return dens_addendum;
}


#endif // NEAL8_IMP_HPP
