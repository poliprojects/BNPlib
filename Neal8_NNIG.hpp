#include <tuple>
#include <vector>
#include <armadillo>
#include <Eigen/Dense> 
#include <stan/math/prim/mat.hpp>
#include <type_traits>
#include "includes.hpp"
#include <math.h>   
// N-NIG model == gaussian kernel + N-IG base measure:
// f ~ N(mu,sig^2)
// (mu,sig^2) ~ G
// G ~ DP(M, G0)  with G0 = N-IG


// Normal likelihoood, Normal Inverse Gamma hierarchy
template<class Hypers> //Hypers = TupleWrapper, distro, ...
class NNIGHierarchy {
protected:
    using state_t= std::array<par_t,2>;

    std::mt19937  rng;
    state_t state; // current values for F's parameters: mu, sigma


    std::shared_ptr<Hypers> hypers; // current values for G0's parameters:mu_0,Lambda0, alpha, beta

public:
    // Contructors:
    ~NNIGHierarchy() = default;
    NNIGHierarchy(std::shared_ptr<Hypers> hypers):
    hypers(hypers)  {}

    // Getters/setters:
    state_t get_state(){return state;}
    void set_state(const state_t &s){state = s;}
    void set_state(int pos, par_t val){state[pos] = val;}

    int get_count(){return hypers.use_count();}



    double log_like(data_t datum) {
        return exp(stan::math::normal_lpdf(datum, state[0], state[1]));
    }

    void draw() {
        float sigmaNew = stan::math::inv_gamma_rng(hypers->get_alpha0(), hypers->get_beta0(), rng);
        float muNew = stan::math::normal_rng(hypers->get_mu0(), sigmaNew/hypers->get_lambda(), rng);
        state[0] = muNew;
        state[1] = sigmaNew;
        }

    void sample_given_data(std::vector<data_t> data) {
        // Get current values of parameters
        auto mu0    = hypers->get_mu0();
        auto Lambda0  = hypers->get_lambda();
        auto alpha0 = hypers->get_alpha0();
        auto beta0  = hypers->get_beta0();

        arma::vec temp = normalGammaUpdate(
          data, mu0, alpha0, beta0, Lambda0);


        auto mu_post = temp(0);
        auto alpha_post = temp(1);
        auto beta_post = temp(2);
        auto postLambda = temp(3);

        // Get a sample
        par_t sigma_new = stan::math::inv_gamma_rng(alpha_post, beta_post, rng);
        par_t mu_new = stan::math::normal_rng(mu_post, sigma_new/postLambda, rng); //? is it ok /postLambda?
        state[0] = mu_new;
        state[1] = sigma_new;
    }



  arma::vec normalGammaUpdate(arma::vec data, double mu0, double alpha0, double beta0,double Lambda0) {
    	double mu_post, alpha_post, beta_post, postLambda;
    	int n = data.size();
    	if (n == 0) {
    		return arma::vec{mu0, alpha0, beta0, Lambda0};
  	}
  	double ybar = arma::mean(data);
  	mu_post = (Lambda0 * mu0 + n * ybar) / (Lambda0 + n);
  	alpha_post = 1.0 * alpha0 + 1.0 * n / 2;

  	// arma::var(x, 1) divides by n, not n-1
  	double ss = n * arma::var(data, 1);

  	beta_post = (beta0 + 0.5 * ss + 0.5 * Lambda0 / (n + Lambda0) * n * std::pow((ybar - mu0), 2));

  	postLambda = Lambda0 + n;

  	return arma::vec{mu_post, alpha_post, beta_post, postLambda};
  }

};


template<template <class> class Hierarchy, class Hypers, class Mixture> // TODO change to MixingMode?
class Neal8{
private:
    unsigned int n_aux=3;
    unsigned int maxiter = 500; // TODO LATER
    unsigned int burnin = 0;
    std::mt19937 rng;
    int numClusters;
    Mixture mixture;
    //Hypers hy;
    //arma::vec probas;


    std::vector<data_t> data;
    std::vector<unsigned int> allocations; // the c vector
    std::vector<Hierarchy<Hypers>> unique_values;
    std::vector<Hierarchy<Hypers>> aux_unique_values;



    void initalize(){
   std::default_random_engine generator;
   std::uniform_int_distribution<int> distribution(0,numClusters);
	
    for (int h = 0; h < numClusters; h++) {
          allocations.push_back(h);
        }
   
      for (int j = numClusters; j < data.size(); j++) {
          int num = distribution(generator); //da stan?
          allocations[j] = num;
        }
    }


    void step(){
        sample_allocations();
        sample_unique_values();
    }

    void sample_allocations(){
        // Other ideas:
        // * our own for loop for k and bool (ci is a singleton)
        // * function from std count distinct values in vector
        // * using a (multi)map?

        // Initialize some relevant variables
        unsigned int k, n_unique, singleton;
	unsigned int n=data.size();
	
        // Initialize cardinalities of unique values

        
 
        for(int i=0; i<n; i++){ // for each data unit data[i]

std::vector<int> card(unique_values.size(), 0);      
for(int j=0; j<n; j++)
            card[ allocations[j] ] += 1;
    
            singleton = 0;
            n_unique = unique_values.size(); //


            if(card[ allocations[i] ] == 1){ // datum i is a singleton
                k = n_unique - 1;
                aux_unique_values[0].set_state( unique_values[allocations[i]
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
            //arma::vec probas(k+n_aux);
	    auto M = mixture.get_totalmass();
	    double tot=0.0;
            for(int k=0; k<n_unique ; k++){ // if datum i is a singleton, then
                // card[k] when k=allocations[i] is equal to 0 -> probas[k]=0

                
                // TODO LATER "meglio in logscale" (?)
                probas(k,0) = card[k] * unique_values[k].log_like(data[i]) / (
                    n-1+M);
		tot+=probas(k,0);
              
            }
            

            for(int k=0; k<n_aux ; k++){
                probas(n_unique+k,0) = (M/n_aux) *
                    aux_unique_values[k].log_like(data[i]) / (n-1+M);
		tot+=probas(n_unique+k,0);
               }
	     probas=probas*(1/tot);

	 std::cout<<"NUNIQUE"<<n_unique<<std::endl;   

            unsigned int c_new = stan::math::categorical_rng(probas, rng);
            if(singleton == 1){
                if (c_new >= n_unique){ // case 1 of 4: SINGLETON - AUX
std::cout<<"test1"<<std::endl;

                    unique_values[ allocations[i] ].set_state(
                        aux_unique_values[c_new-n_unique].get_state());
                    card[ allocations[i] ] += 1;
                }
                else{ // case 2 of 4: SINGLETON - OLD VALUE
std::cout<<"test2"<<std::endl;

                    unique_values.erase(
                        unique_values.begin()+allocations[i] ); // TODO ho levato il -1, giusto? allocations parte da 0

                    card.erase( card.begin()+allocations[i] );
                    card[c_new] += 1;

                    int tmp = allocations[i];
std::cout<<"c_new"<<c_new<<std::endl;
std::cout<<"alloc[i]"<<tmp<<std::endl;

                    allocations[i] = c_new;
                    for(auto &c : allocations){ // relabeling
                        if (c > tmp)
                            c -= 1;
                    }



for(unsigned int i=0; i<n; i++){ 
 std::cout<<"i:"<<i<<std::endl;
 std::cout<<"allo:"<<allocations[i]<<std::endl; }

                }
            } // end of if(singleton == 1)

            else{ // if singleton == 0

                if (c_new>=n_unique){ // case 3 of 4: NOT SINGLETON - AUX
std::cout<<"test3"<<std::endl;
                    unique_values.push_back(aux_unique_values[c_new-n_unique]);
                    card.push_back(1);
                    allocations[i] = n_unique ; // TODO tolto il -1, ho aggiunto un elem

                }
                else{ // case 4 of 4: NOT SINGLETON - OLD VALUES
std::cout<<"test4"<<std::endl;

                    allocations[i] = c_new;
                    card[c_new] += 1;
                }

            }
        } // end of for(int i=0; i<n; i++) loop

    } // end of sample_allocations()



    void sample_unique_values(){

        numClusters=unique_values.size();

        std::vector<std::vector<unsigned int>> clust_idxs(numClusters);
	unsigned int n=allocations.size();


        for(unsigned int i=0; i<n; i++){ // save different cluster in each row
             clust_idxs[ allocations[i]].push_back(i);}



	 for(int j=0; j<numClusters; j++){ 
std::cout<<"cluster num "<< j;
             for (unsigned int i=0; i<clust_idxs[j].size(); i++ )
            std::cout<<" "<<clust_idxs[j][i];
std::cout<<""<<std::endl;
}

        for (unsigned int j=0; j< numClusters; j++) {
            std::vector<data_t> curr_data;

            for ( auto &idx : clust_idxs[j] )
                curr_data.push_back( data[idx] );

            unique_values[j].sample_given_data(curr_data);
        }


    std::cout<<""<<std::endl;
    }



    void save_iteration(unsigned int iter){ // TODO LATER
        std::cout << "Iteration n. " << iter << " / " << maxiter << std::endl;
        print();
    }

    void print() {
        for (int h = 0; h < numClusters; h++) {
            std::cout << "Cluster # " << h << std::endl;
	    for (auto c:unique_values[h].get_state()){
	            std::cout << "Parameters: "<< c<<std::endl;
            }
          }
        }



public:

    ~Neal8() = default;
    Neal8(const std::vector<data_t> & data, int numClusters,int n_aux, const Mixture & mix,const Hypers &hy ):
    data(data), numClusters(numClusters), n_aux(n_aux), mixture(mix)  {
	
      Hierarchy<Hypers> hierarchy(std::make_shared<Hypers> (hy));
      for (int h = 0; h < numClusters; h++) {
              unique_values.push_back(hierarchy);
            }
      for (int h = 0; h < n_aux; h++) {
        aux_unique_values.push_back(hierarchy);
      }
	//std::cout<<unique_values[0].get_count();
    }


    Neal8(std::vector<data_t>  & data, int n_aux, const Mixture & mix, const Hypers &hy):
    Neal8(data, data.size(), n_aux, mix, hy ) {}


    void run(){
        initalize();
 		unsigned int iter = 0;
        while(iter < maxiter){
            step();    
            if(iter >= burnin)
              save_iteration(iter);
			iter++;
        }
    }

}; // end of Class Neal8


//TODO LIST and general doubts

// - is it correct the posterior update given data, given clusters?

// - a way to erase and add clusters more efficiently

class SimpleMixture {
double totalmass;

public:

    ~SimpleMixture() = default;
    SimpleMixture(double totalmass):
    totalmass(totalmass) {
      assert(totalmass>=0);
    }

    double const get_totalmass(){return totalmass;}
};


class HypersFixed {
    double mu0, lambda, alpha0, beta0;
public:
    double get_mu0(){return mu0;}
    double get_alpha0(){return alpha0;}
    double get_beta0(){return beta0;}
    double get_lambda(){return lambda;}

    ~HypersFixed() = default;
    HypersFixed(double mu0, double lambda, double alpha0, double beta0):
    mu0(mu0), lambda(lambda), alpha0(alpha0), beta0(beta0)  {}
};
