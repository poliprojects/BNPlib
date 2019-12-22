// DISTRIBUTION WRAPPER
template <class DistributionType, class... HyperDistributionTypes>
class DistributionWrapper{
    private:
        //arma::vec parameters;
        DistributionType distribution;
        bool has_prior;
        std::tuple<HyperDistributionTypes...> priors;
    public:
        DistributionWrapper(const std::tuple<HyperDistributionTypes...> &args) :
            priors(args) {}
        // Other constructors

        //data_t sample_from_priors() {...}
        //data_t sample(){
        //  ...
        //  std::vector<data_t> samples;
        //  for(int i = 0; i < N; i++){
        //      val = priors[i].sample();
        //  }
        //}
};

// Vecchia idea:
// class DistributionWrapper{
//  private:
//      arma::vec parameters;
//      shared_ptr<DistributionWrapper> prior = nullptr;
//      arma::vec density() virtual;
//  public:
//      arma::vec sample_from_prior();
//      Distribution(arma::vec parameters)
// }

// Idea by Mario:
template<state_t, data_t>
class DistributionWrapper{

    real eval(data_t x){
        ...
    }
    void update(std::vector<state_t> vals){
        ...
    }
}




// HIERARCHY
// Lik e P0 ereditano da DistributionWrapper
template <typename Lik, typename P0>
class Hierarchy {
    Lik lik;
    P0* po;

    public:
        typename data_t
        sample(std::vector<data_t> data); // specifico della hierarchy

        real log_lik(data_t datum) {
            return lik.eval(datum);
        }
}

class HypersFixed {
    double mu0, k0, alpha0, beta0;

    //getters

    void sample() {return;}

}

class Hypers{
    real current_value;
    void sample(){...} // prende a caso
}









// DISTRIBUTION
class Distribution{
	arma::vec parameters;
	double eval() virtual;
}

class GaussianDistribution: public Distribution{
    double eval() override{
			return 1/(2*pi)^0.5 * ...;
		}
}
// Poi creo oggetti gaussian, etc:
// Distribution gaussian;


class DistributionWithHypers{
	std::vector<shared_ptr<DistributionWithFixed>> hyperparameters;
	arma::vec sample();
}

class DistributionWithFixed: public DistributionWithHypers{
	arma::vec parameters;




// HIERARCHY
class Hierarchy{
	DistributionWithHypers base_measure;
	DistributionWithHypers kernel;
}
