#ifndef NEAL2_NNIG_IMP_HPP
#define NEAL2_NNIG_IMP_HPP

#include "Neal2_NNIG.hpp"

// TODO transfer over all changes from Neal8 to here

template<template <class> class Hierarchy, class Hypers, class Mixture>
void Neal2<Hierarchy,Hypers,Mixture>::initalize(){
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0,num_clusters);

    for (int h = 0; h < num_clusters; h++) {
      allocations.push_back(h);
    }
    for (int j = num_clusters; j < data.size(); j++) {
        int num = distribution(generator); //TODO da stan?
        allocations[j] = num;
    }
    }

template<template <class> class Hierarchy, class Hypers, class Mixture>
void Neal2<Hierarchy,Hypers,Mixture>::sample_allocations(){
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
        for(int j = 0; j < n; j++)
            card[ allocations[j] ] += 1;

        singleton = 0;
        n_unique = unique_values.size();

        if(card[ allocations[i] ] == 1){ // datum i is a singleton
            singleton = 1;
        }
        
        card[ allocations[i] ] -= 1;

        
        // Draw a NEW value for ci
        Eigen::MatrixXd probas(n_unique +(1-singleton) ,1); 
        

        auto M = mixture.get_totalmass();
        double tot=0.0;

        for(int k = 0; k < n_unique ; k++){
		
            probas(k,0) = card[k] * unique_values[k].log_like(data[i]) / (
                n-1+M);
			if(singleton==1 && k==i){
				// Take the hyperparameters
				std::shared_ptr<Hypers> hy= unique_values[0].get_hypers();
				double mu0    = hy->get_mu0();
				double lambda = hy->get_lambda();
				double alpha0 = hy->get_alpha0();
				double beta0  = hy->get_beta0();

				// Compute pieces of the weight
				//double factor1 = alpha0 * sqrt(lambda) * pow(2*beta0,alpha0) *
					//tgamma(alpha0/2 + 0.25);
				//double factor2 = sqrt(M_PI) * pow(1+lambda, alpha0+1) *
					//tgamma(alpha0+1);
				//double y = 1.0; // TODO what is y?
				//double base = (y*y + lambda*mu0*mu0+2*beta0)/(1+lambda) -
					//(y+mu0)*(y+mu0) / ((1+lambda)*(1+lambda));
				
				// Compute the weight
				//probas(i,0) = factor1 / factor2 * pow(base, 1.5);
				double sigtilde= sqrt(beta0*(lambda+1/(alpha0*lambda)));
				probas(i,0) = M * exp(stan::math::student_t_lpdf(data[i],
                    2*alpha0, mu0, sigtilde))/ (n-1+M);

			} // metti in i la probas c!=cj con integrale TODO done?
            tot+=probas(k,0);
        }

		if(singleton == 0){
			std::shared_ptr<Hypers> hy= unique_values[0].get_hypers();
				double mu0    = hy->get_mu0();
				double lambda = hy->get_lambda();
				double alpha0 = hy->get_alpha0();
				double beta0  = hy->get_beta0();

			double sigtilde= sqrt(beta0*(lambda+1/(alpha0*lambda)));
			probas(n_unique,0) = M * exp(stan::math::student_t_lpdf(data[i],
                2*alpha0, mu0, sigtilde))/ (n-1+M);
			tot+=probas(n_unique,0);
		}

        
        probas = probas * (1/tot);

  		
        
        unsigned int c_new = stan::math::categorical_rng(probas, rng) - 1;
        
       

        if(singleton == 1){
            if(c_new == allocations[i]){ // case 1 of 4: SINGLETON - SINGLETON
            	std::shared_ptr<Hypers> hy= unique_values[0].get_hypers();
				double mu0    = hy->get_mu0();
				double lambda = hy->get_lambda();
				double alpha0 = hy->get_alpha0();
				double beta0  = hy->get_beta0();

				std::array<double,2> par_pair;
				double sigma_new = stan::math::inv_gamma_rng(alpha0 + 1/2,
        			beta0+(lambda*pow(data[i]-mu0,2))/(lambda+2), rng);
    			double mu_new = stan::math::normal_rng((lambda*mu0+data[i])/(lambda+1),
        			sigma_new/(lambda+1), rng); 
				par_pair[0]=mu_new;
				par_pair[1]=sigma_new;
				unique_values[ allocations[i]].set_state(par_pair); // TODO draw da H
				
            }
            else{ // case 2 of 4: SINGLETON - CLUSTER
                unique_values.erase(
                    unique_values.begin()+allocations[i] );
                
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
            if(c_new == n_unique){ // case 3 of 4: NOT SINGLETON - SINGLETON
				std::shared_ptr<Hypers> hy= unique_values[0].get_hypers();
				double mu0    = hy->get_mu0();
				double lambda = hy->get_lambda();
				double alpha0 = hy->get_alpha0();
				double beta0  = hy->get_beta0();

				std::array<double,2> par_pair;
				double sigma_new = stan::math::inv_gamma_rng(alpha0 + 1/2,
        			beta0+(lambda*pow(data[i]-mu0,2))/(lambda+2), rng);
    			double mu_new = stan::math::normal_rng((lambda*mu0+data[i])/(lambda+1),
        			sigma_new/(lambda+1), rng); 
				par_pair[0]=mu_new;
				par_pair[1]=sigma_new;
				Hierarchy<Hypers> newUnique(hy);
				newUnique.set_state(par_pair); 
                unique_values.push_back(newUnique); 
                allocations[i] = n_unique;
            }
            else{ // case 4 of 4: NOT SINGLETON - CLUSTER
                allocations[i] = c_new;
            }
        } // end of else

    } // end of for(int i = 0; i < n; i++) loop

    } // end of sample_allocations()




template<template <class> class Hierarchy, class Hypers, class Mixture>
void Neal2<Hierarchy,Hypers,Mixture>::sample_unique_values(){

    num_clusters=unique_values.size();
    std::vector<std::vector<unsigned int>> clust_idxs(num_clusters);
    unsigned int n = allocations.size();
    for(int i = 0; i < n; i++){ // save different cluster in each row
        clust_idxs[ allocations[i] ].push_back(i);
    }

    // DEBUG:
    for(int j = 0; j < num_clusters; j++){ 
        std::cout << "Cluster #" << j << ": ";
        for (unsigned int i = 0; i < clust_idxs[j].size(); i++)
            std::cout << " " << clust_idxs[j][i];
        std::cout << std::endl;
    }

    for (int j = 0; j < num_clusters; j++) {
        std::vector<double> curr_data;
        for ( auto &idx : clust_idxs[j] )
            curr_data.push_back( data[idx] );
        unique_values[j].sample_given_data(curr_data);
    }

    std::cout << std::endl; // DEBUG
    }


template<template <class> class Hierarchy, class Hypers, class Mixture>
void Neal2<Hierarchy,Hypers,Mixture>::save_iteration(unsigned int iter){
    // TODO
    std::cout << "Iteration n. " << iter+1 << " / " << maxiter << std::endl;
    print();
    }


template<template <class> class Hierarchy, class Hypers, class Mixture>
void Neal2<Hierarchy,Hypers,Mixture>::print(){
    for (int h = 0; h < num_clusters; h++) {
        std::cout << "Cluster # " << h << std::endl;
        std::cout << "Parameters: ";

        for (auto c : unique_values[h].get_state()){
            std::cout << c << " " << std::endl;
        }
        std::cout << std::endl;
    }
    }




#endif // NEAL2NNIG_HPP
