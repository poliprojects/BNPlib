#ifndef NEAL8_IMP_HPP
#define NEAL8_IMP_HPP

#include "Neal8.hpp"

template<template <class> class Hierarchy, class Hypers, class Mixture>
const void Neal8<Hierarchy, Hypers, Mixture>::print_startup_message(){
    std::cout << "Running Neal8 algorithm (with m=" << n_aux <<
        " auxiliary blocks)..." << std::endl;
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
void Neal8<Hierarchy, Hypers, Mixture>::sample_allocations(){
    
    // TODO Other ideas:
    // * our own for loop for k and bool (ci is a singleton)
    // * function from std count distinct values in vector
    // * using a (multi)map?
    // Initialize some relevant variables
    unsigned int k, n_unique, singleton;
    unsigned int n = this->data.rows();
  
    for(int i = 0; i < n; i++){ // for each data unit data[i]
        singleton = 0;
        n_unique = this->unique_values.size();
     
        if(this->cardinalities[ this->allocations[i] ] == 1){
            // datum i is a singleton
            k = n_unique - 1; // TODO ci serve ?
            aux_unique_values[0].set_state( this->unique_values[
                this->allocations[i] ].get_state() ); // move phi value in aux
            singleton = 1;
        }
        else{
            k = n_unique;
        }

        // Remove point from cluster
        this->cardinalities[ this->allocations[i] ] -= 1;
        
        // Draw the aux from G0
        for(int j = singleton; j < n_aux; j++){
            aux_unique_values[j].draw();
        }

        // Compute probabilities of extracting each cluster
        Eigen::VectorXd probas(n_unique+n_aux); //k or n_unique (TODO ?)
        
        double tot = 0.0;
        for(int k = 0; k < n_unique; k++){ // if datum i is a singleton, then
            // card[k] when k=allocations[i] is equal to 0 -> probas[k]=0
            probas(k) = this->mixture.prob_existing_cluster(
                this->cardinalities[k], n ) *
                this->unique_values[k].like(this->data.row(i))(0);
            tot += probas(k);
        }

        for(int k = 0; k < n_aux; k++){
            probas(n_unique+k) = this->mixture.prob_new_cluster(n, n_unique) *
                aux_unique_values[k].like(this->data.row(i))(0) / n_aux;
            tot += probas(n_unique+k);
        }

        // Normalize
        probas = probas / tot;

        // Draw a NEW value for ci
        unsigned int c_new = stan::math::categorical_rng(probas, this->rng) - 1;

        if(singleton == 1){
            if(c_new >= n_unique){ // case 1 of 4: SINGLETON - AUX
                this->unique_values[ this->allocations[i] ].set_state(
                    aux_unique_values[c_new-n_unique].get_state());
                this->cardinalities[ this->allocations[i] ] += 1;
            }
            else{ // case 2 of 4: SINGLETON - OLD VALUE
                this->unique_values.erase(
                    this->unique_values.begin() + this->allocations[i] );

                this->cardinalities.erase(
                	this->cardinalities.begin()+this->allocations[i] );
                this->cardinalities[c_new] += 1;

                unsigned int c_old = this->allocations[i];

                this->allocations[i] = c_new;
                for(auto &c : this->allocations){ // relabeling
                    if(c > c_old){
                        c -= 1;
                    }
                }
            } // end of else
        } // end of if(singleton == 1)
        
        else{ // if singleton == 0
            if(c_new >= n_unique){ // case 3 of 4: NOT SINGLETON - AUX
                this->unique_values.push_back(
                    aux_unique_values[c_new-n_unique]);
                this->allocations[i] = n_unique;
                this->cardinalities.push_back(1);
            }
            else{ // case 4 of 4: NOT SINGLETON - OLD VALUES
                this->allocations[i] = c_new;
                this->cardinalities[c_new] += 1;
            }
        } // end of else

    } // end of for(int i = 0; i < n; i++) loop

} // end of sample_allocations()


template<template <class> class Hierarchy, class Hypers, class Mixture>
Eigen::VectorXd Neal8<Hierarchy, Hypers, Mixture>::density_marginal_component(
    Hierarchy<Hypers> &temp_hier, unsigned int n){
    double M = this->mixture.get_totalmass();
    Eigen::VectorXd dens_addendum(this->density.first.rows());
            
    // Component from G0
    for(int h = 0; h < n_aux; h++){
        temp_hier.draw();
        dens_addendum += temp_hier.like(this->density.first) * (M/n_aux)/(M+n);
    }

    return dens_addendum;
}


#endif // NEAL8_IMP_HPP
