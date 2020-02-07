#ifndef NEAL8_NNIG_IMP_HPP
#define NEAL8_NNIG_IMP_HPP

#include <tuple>
#include <vector>
#include <Eigen/Dense>
#include <stan/math/prim/mat.hpp>
#include <type_traits>
#include "includes_universal.hpp"
#include <math.h>
#include "NNIGHierarchy.hpp"
#include "SimpleMixture.hpp"
#include "HypersFixed.hpp"
// N-NIG model == gaussian kernel + N-IG base measure:
// f ~ N(mu,sig^2)
// (mu,sig^2) ~ G
// G ~ DP(M, G0)  with G0 = N-IG
// #include "output.pb.h"
#include "Neal8_NNIG.hpp"
// Normal likelihoood, Normal Inverse Gamma hierarchy


template<template <class> class Hierarchy, class Hypers, class Mixture>
void Neal8<Hierarchy,Hypers,Mixture>::initalize(){
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0,numClusters);

    for (int h = 0; h < numClusters; h++) {
      allocations.push_back(h);
    }
    for (int j = numClusters; j < data.size(); j++) {
        int num = distribution(generator); //TODO da stan?
        allocations[j] = num;
    }
    }

template<template <class> class Hierarchy, class Hypers, class Mixture>
void Neal8<Hierarchy,Hypers,Mixture>::sample_allocations(){
    // TODO Other ideas:
    // * our own for loop for k and bool (ci is a singleton)
    // * function from std count distinct values in vector
    // * using a (multi)map?
    // Initialize some relevant variables
    unsigned int k, n_unique, singleton;
    unsigned int n=data.size();

    for(int i=0; i<n; i++){ // for each data unit data[i]

        // Initialize cardinalities of unique values
        std::vector<int> card(unique_values.size(), 0);
        for(int j=0; j<n; j++)
            card[ allocations[j] ] += 1;

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
        for(int j=singleton; j<n_aux; j++){
            aux_unique_values[j].draw();
        }

        // Draw a NEW value for ci
        Eigen::MatrixXd probas(n_unique+n_aux,1); //k or n_unique
        //Matrix<double, Dynamic, 1> VectorXd

        auto M = mixture.get_totalmass();
        double tot=0.0;
        for(int k=0; k<n_unique ; k++){ // if datum i is a singleton, then
            // card[k] when k=allocations[i] is equal to 0 -> probas[k]=0

            // TODO LATER "meglio in logscale" (?)
            probas(k,0) = card[k] * unique_values[k].log_like(data[i]) / (
                n-1+M);
            tot+=probas(k,0);
        }

        for(int k=0; k<n_aux; k++){
            probas(n_unique+k,0) = (M/n_aux) *
                aux_unique_values[k].log_like(data[i]) / (n-1+M);
            tot += probas(n_unique+k,0);
           }
        probas = probas * (1/tot);

        //for(int i=0; i<probas.size(); i++){
        //    std::cout << "probas_" << probas(i,0) << std::endl; // DEBUG
        //}
        unsigned int c_new = stan::math::categorical_rng(probas, rng) -1;

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
                    if (c > tmp)
                        c -= 1;
                }
            } // end of else
        } // end of if(singleton == 1)
        else{ // if singleton == 0
            if (c_new>=n_unique){ // case 3 of 4: NOT SINGLETON - AUX
                unique_values.push_back(aux_unique_values[c_new-n_unique]);
                card.push_back(1);
                allocations[i] = n_unique;
            }
            else{ // case 4 of 4: NOT SINGLETON - OLD VALUES
                allocations[i] = c_new;
                card[c_new] += 1;
            }
        } // end of else

    } // end of for(int i=0; i<n; i++) loop

    } // end of sample_allocations()




template<template <class> class Hierarchy, class Hypers, class Mixture>
void Neal8<Hierarchy,Hypers,Mixture>::sample_unique_values(){

    numClusters=unique_values.size();
    std::vector<std::vector<unsigned int>> clust_idxs(numClusters);
    unsigned int n = allocations.size();
    for(unsigned int i=0; i<n; i++){ // save different cluster in each row
        clust_idxs[ allocations[i] ].push_back(i);
    }

    // DEBUG:
    //for(int j=0; j<numClusters; j++){
    //    std::cout << "Cluster #" << j << ": ";
    //    for (unsigned int i=0; i<clust_idxs[j].size(); i++)
    //        std::cout << " " << clust_idxs[j][i];
    //    std::cout << std::endl;
    //}

    for (unsigned int j=0; j< numClusters; j++) {
        std::vector<data_t> curr_data;
        for ( auto &idx : clust_idxs[j] )
            curr_data.push_back( data[idx] );
        unique_values[j].sample_given_data(curr_data);
    }

    // std::cout << std::endl; // DEBUG
    }


template<template <class> class Hierarchy, class Hypers, class Mixture>
void Neal8<Hierarchy,Hypers,Mixture>::save_iteration(unsigned int iter){
	// TODO

    //DEBUG:
    //std::cout << "Iteration # " << iter << " / " << maxiter-1 << std::endl;
	IterationOutput iter_out;

    *iter_out.mutable_allocations() = {allocations.begin(), allocations.end()};

	for (int i=0; i<unique_values.size(); i++){
		UniqueValues temp;
        for(auto &par : unique_values[i].get_state())
            temp.add_params(par);
        iter_out.add_phi();
        *iter_out.mutable_phi(i) = temp;
	}

	chain.add_state();
	*chain.mutable_state(iter) = iter_out;

    //print();
}

template<template <class> class Hierarchy, class Hypers, class Mixture>
void Neal8<Hierarchy,Hypers,Mixture>::cluster_estimate(){
	
	Eigen::VectorXf errors(maxiter);


	int n=data.size();
	
	Eigen::MatrixXd TOTdiss(n,n);
	TOTdiss=Eigen::MatrixXd::Zero(n, n);
	std::vector<Eigen::MatrixXd> ALLdiss;
	IterationOutput temp;
	
	for(int h = 0; h < maxiter; h++){
		temp= *chain.mutable_state(h);
		Eigen::MatrixXd dissim(n,n);
		dissim=Eigen::MatrixXd::Zero(n, n);
		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				if(temp.allocations(i)==temp.allocations(j))
				dissim(i,j) =1;
			}
		}

		ALLdiss.push_back(dissim);
		TOTdiss=TOTdiss+dissim;
	
	}

	TOTdiss=TOTdiss/maxiter;


	for(int h = 0; h < maxiter; h++){

		errors(h)=(TOTdiss-ALLdiss[h]).norm();
		
	}
	
	//std::cout<<errors<<std::endl; //DEBUG
	std::ptrdiff_t i;
	int minerr = errors.minCoeff(&i);
	//std::cout<<i << " "; //DEBUG
	//return chain[i];
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
void Neal8<Hierarchy,Hypers,Mixture>::eval_density(){
Eigen::VectorXf density(grid.size());
auto M = mixture.get_totalmass();
int n=data.size();
IterationOutput temp;
for (int i = 0; i < maxiter; i++) {
	temp= *chain.mutable_state(i);
	std::vector<int> card(temp.phi_size(), 0); // TODO salviamoci ste card da qualche parte
	for(int j=0; j<n; j++)
        card[ temp.allocations(j) ] += 1;

	Hierarchy<Hypers> temp_uniq_v(unique_values[0].get_hypers());
	for(int h = 0; h < temp.phi_size(); h++) { 
		std::array<par_t,2> temp_state;
		temp_state[0]=temp.phi(h).params(0);
		temp_state[1]=temp.phi(h).params(1);
		temp_uniq_v.set_state(temp_state);
		density+=card[h]*temp_uniq_v.log_like(grid)/(data.size()+M);
	}
	density+=M*temp_uniq_v.log_like(grid)/(M+data.size());
}

density=density/maxiter;
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
void Neal8<Hierarchy,Hypers,Mixture>::print(){
    for (int h = 0; h < numClusters; h++) {
        std::cout << "Cluster # " << h << std::endl;
        std::cout << "Parameters: ";

        for (auto c : unique_values[h].get_state()){
            std::cout << c << " " << std::endl;
        }
        std::cout << std::endl;
    }
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
void Neal8<Hierarchy,Hypers,Mixture>::write_clusters_to_file(
    std::string filename){
    std::ofstream file;
    file.open(filename);

    file << "number,datum,cluster,mu,sigma2" << std::endl;
    for(int i=0; i<data.size(); i++){
        auto params = unique_values[ allocations[i] ].get_state();
        file << i << "," << data[i] << "," << allocations[i] << "," <<
        params[0] << "," << params[1] << std::endl;
    }
    
    file.close();
}



#endif // NEAL8NNIG_HPP
