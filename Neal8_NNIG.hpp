// N-NIG model == gaussian kernel + N-IG base measure:
// f ~ N(mu,sig^2)
// (mu,sig^2) ~ G
// G ~ DP(M, G0)  with G0 ~ N-IG


// Normal likelihoood, Normal Inverse Gamma hierarchy
template<class Hypers> //Hypers = TupleWrapper, distro, ...
class NNIGHierarchy {
protected:
    using std::tuple<par_t, par_t> = state_tuple_t;

    unsigned int rng = 20191225;
    state_tuple_t state; // current values for F's parameters: mu, sigma


    std::shared_ptr<Hypers> hypers; // current values for G0's parameters:mu_0,Lambda0, alpha, beta

public:
    // Contructors:
    // TODO!

    // Getters/setters:
    state_tuple_t get_state(){return state;}
    void set_state(const state_tuple_t &s){state = s;}
    void set_state(int pos, par_t val){state(pos) = val;}

    double log_like(data_t datum) {
        return stan::math::normal_lpdf(datum, state(0), state(1));
    }

    void draw() {
        real sigmaNew = stan::math::inv_gamma_rng(hypers.get_alpha0(),
            hypers.get_beta0(), rng);
        real muNew = stan::math::normal_rng(hypers.get_m0(), sigmaNew/hypers.get_lambda(),
            rng);
        state(0) = muNew;
        state(1) = sigmaNew;
        }

    void sample_given_data(std::vector<data_t> data) {
        // Get current values of parameters
        auto mu0    = hypers.get_mu0();
        auto Lambda0  = hypers.get_lambda();
        auto alpha0 = hypers.get_alpha0();
        auto beta0  = hypers.get_beta0();

        arma::vec temp = normalGammaUpdate(
          data, mu0, alpha0, beta0, Lambda0);

        auto mu_post = temp(0);
        auto alpha_post = temp(1);
        auto beta_post = temp(2);
        auto postLambda = temp(3);

        // Get a sample
        par_t sigma_new = stan::math::inv_gamma_rng(alpha_post, beta_post, rng);
        par_t mu_new = stan::math::normal_rng(mu_post, sigma_new/postLambda, rng); //? is it ok /postLambda?
        state(0) = mu_new;
        state(1) = sigma_new;
    }

}

arma::vec NNIGHierarchy::normalGammaUpdate(
    arma::vec data, double mu0, double alpha0, double beta0,
    double Lambda0) {
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

  beta_post = (
      beta0 + 0.5 * ss +
      0.5 * Lambda0 / (n + Lambda0) * n * std::pow((ybar - mu0), 2));

  postLambda = Lambda0 + n;

  return arma::vec{mu_post, alpha_post, beta_post, postLambda};
}



template<class Hierarchy<Hypers>, unsigned int n_aux = 3>
class Neal8{
private:
    unsigned int m = n_aux; // TODO
    unsigned int maxiter = 1000;
    unsigned int burnin = 0;
    int numClusters;
    //arma::vec probas;


    //std::vector<Hierarchy> hierarchies;
    std::vector<data_t> data;
    std::array<n, unsigned int> allocations; // the c vector
    std::vector<Hierarchy> unique_values;
    std::array<m, Hierarchy> aux_unique_values;

    //Hypers hypers;


    void initalize(){

    for (int h = 0; h < numClusters; h++) {
          Hierarchy hierarchy;
          unique_values.push_back((hierarchy));
          allocations[h] = h;
        }

        for (int j = numClusters; j < data.size(); j++) {
          int num = stats::rdiscreteunif(0, numClusters, engine); //da stan?
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

        // Initialize cardinalities of unique values
        std::vector<unsigned int> card(unique_values.size(), 0);
        for(int j=0; j<n; j++)
            card[ allocations[j] ] += 1;

        for(int i=0; i<n; i++){ // for each data unit data[i]

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
            for(int j=singleton; j<m; j++){
                aux_unique_values[j].draw();
            }

            // Draw a NEW value for ci
            Eigen::Matrix<double, k+m, 1> probas;

            for(int k=0; k<n_unique ; k++){ // if datum i is a singleton, then
                // card[k] when k=allocations[i] is equal to 0 -> probas[k]=0

                auto M = mixture.get_totalmass();

                // TODO giusto? "meglio in logscale" (?)
                probas(k) = card[k] * unique_values[k].loglike(data[i]) / (
                    n-1+M);
            }
            for(int k=0; k<m ; k++){
                probas(n_unique+k) = (M/m) *
                    aux_unique_values[k].loglike(data[i]) / (n-1+M);
            }

            unsigned int c_new = stan::math::categorical_rng(probas,rng);
            if(singleton == 1){
                if (c_new >= n_unique){ // case 1 of 4: SINGLETON - AUX
                    unique_values[ allocations[i] ].set_state(
                        aux_unique_values[c_new-n_unique].get_state());
                    card[ allocations[i] ] += 1;
                }
                else{ // case 2 of 4: SINGLETON - OLD VALUE
                    unique_values.erase(
                        unique_values.begin()+allocations[i]-1 );

                    card.erase( card.begin()+allocations[i]-1 );
                    card[c_new] += 1;

                    int tmp = allocations[i];
                    allocations[i] = c_new;
                    for(auto &c : allocations){ // relabeling
                        if (c > tmp)
                            c -= 1;
                    }
                }
            } // end of if(singleton == 1)

            else{ // if singleton == 0

                if (c_new>=n_unique){ // case 3 of 4: NOT SINGLETON - AUX
                    unique_values.push_back(aux_unique_values[c_new-n_unique]);
                    card.push_back(1);
                    allocations[i] = n_unique + 1;
                }
                else{ // case 4 of 4: NOT SINGLETON - OLD VALUES
                    allocations[i] = c_new;
                    card[c_new] += 1;
                }

            }
        } // end of for(int i=0; i<n; i++) loop

    } // end of sample_allocations()



    void sample_unique_values(){
        numClusters=unique_values.size();

        std::vector<std::vector<unsigned int>> clust_idxs;
        for(int i=0; i<n; i++) // save different cluster in each row
            clust_idxs[ allocations[i] ].push_back(i);

        for (int j=0; j< clust_idxs.size(); j++) {
            std::vector<data_t> curr_data;

            for ( auto &idx : clust_idxs[j] )
                curr_data.push_back( data[idx] );

            unique_values[j].sample(curr_data);
        }
    }



    void save_iteration(unsigned int iter){ // TODO
        std::cout << "Iteration n. " << iter << " / " << maxiter << std::endl;
        print();
    }

    void print() {
        for (int h = 0; h < numClusters; h++) {
            std::cout << "Cluster # " << h << std::endl;
            std::cout << "Parameters: "<< unique_values[h].getstate();
          }
        }



public:

    ~Neal8() = default;
    Neal8(std::vector<data_t> & data, int numClusters):
    data(data), numClusters(numClusters) {};
    Neal8(std::vector<data_t>  & data, int m):
    data(data), numClusters(data.size()) {};
    // TODO!

    void run(){
        initalize();
        unsigned int iter = 0;
        while(iter < maxiter){
            step();
            if(iter >= burnin)
                save_iteration(iter);
        }
    }

} // end of Class Neal8


//TO DO LIST and general doubts

// - is it correct the posterior update given data, given clusters?
// - Hierarchy constructors
// - Hypers constructors
// unif con stan
// - a way to erase and add clusters more efficient




class HypersFixed {
    double mu0, lambda, alpha0, beta0;
    double get_mu0(){return mu0;}
    double get_alpha0(){return alpha0;}
    double get_beta0(){return beta0;}
    double get_lambda(){return lambda;}



}


int main() {
    Neal8<NNIGHierarchy<HypersFixed>> sampler;
}
