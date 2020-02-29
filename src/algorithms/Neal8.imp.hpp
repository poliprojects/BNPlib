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
    unsigned int n = data.size();

    for(int i = 0; i < n; i++){ // for each data unit data[i]

        // Initialize cardinalities of unique values
        std::vector<int> card(unique_values.size(), 0);
        for(int j = 0; j < n; j++){
            card[ allocations[j] ] += 1;
        }

        singleton = 0;
        n_unique = unique_values.size();

        if(card[ allocations[i] ] == 1){ // datum i is a singleton
            k = n_unique - 1;
            aux_unique_values[0].set_state( unique_values[ allocations[i]
                ].get_state() ); // move phi value in aux
            singleton = 1;
        }
        else{
            k = n_unique;
        }

        card[ allocations[i] ] -= 1;

        // Draw the aux from G0
        for(int j = singleton; j < n_aux; j++){
            aux_unique_values[j].draw();
        }

        // Draw a NEW value for ci
        Eigen::VectorXd probas(n_unique+n_aux); //k or n_unique

        //auto M = mixture.get_totalmass();
        double tot = 0.0;
        for(int k = 0; k < n_unique ; k++){ // if datum i is a singleton, then
            // card[k] when k=allocations[i] is equal to 0 -> probas[k]=0
            probas(k) = mixture.prob_existing_cluster(card[k],n) *
            	unique_values[k].like(data[i]);
            tot += probas(k);
        }

        for(int k = 0; k < n_aux; k++){
            probas(n_unique+k) = mixture.prob_new_cluster(n, n_unique) *
                aux_unique_values[k].like(data[i]) / n_aux;
            tot += probas(n_unique+k,0);
        }
        probas = probas / tot;

        //for(int i = 0; i < probas.size(); i++){
        //    std::cout << "probas_" << probas(i) << std::endl; // DEBUG
        //}
        unsigned int c_new = stan::math::categorical_rng(probas, rng) - 1;

        //std::cout<<"c_new: "<<c_new<<std::endl; // DEBUG

        if(singleton == 1){
            if(c_new >= n_unique){ // case 1 of 4: SINGLETON - AUX
                unique_values[ allocations[i] ].set_state(
                    aux_unique_values[c_new-n_unique].get_state());
                card[ allocations[i] ] += 1;
            }
            else{ // case 2 of 4: SINGLETON - OLD VALUE
                unique_values.erase(
                    unique_values.begin()+allocations[i] );
                card.erase( card.begin()+allocations[i] );
                card[c_new] += 1;
                int tmp = allocations[i];
                allocations[i] = c_new;
                for(auto &c : allocations){ // relabeling
                    if(c > tmp){
                        c -= 1;
                    }
                }
            } // end of else
        } // end of if(singleton == 1)
        else{ // if singleton == 0
            if(c_new >= n_unique){ // case 3 of 4: NOT SINGLETON - AUX
                unique_values.push_back(aux_unique_values[c_new-n_unique]);
                card.push_back(1);
                allocations[i] = n_unique;
            }
            else{ // case 4 of 4: NOT SINGLETON - OLD VALUES
                allocations[i] = c_new;
                card[c_new] += 1;
            }
        } // end of else

    } // end of for(int i = 0; i < n; i++) loop

} // end of sample_allocations()


template<template <class> class Hierarchy, class Hypers, class Mixture>
unsigned int Neal8<Hierarchy, Hypers, Mixture>::cluster_estimate(){
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


#endif // NEAL8_IMP_HPP
