namespace bnplib{
	using data_t = double;
	using Gamma = std::gamma_distribution<data_t>;
	using Normal = std::normal_distribution<data_t>;
	// etc...
}

using namespace bnplib;


//template <size_t N>
//class unpack_caller{
//	private:
//	    template <typename FuncType, size_t... I>
//	    void call(FuncType &f, std::vector<data_t> &args, indices<I...>){
//	        f(args[I]...);
//	    }
//	public:
//	    template <typename FuncType>
//	    void operator () (FuncType &f, std::vector<data_t> &args){
//	        call(f, args, BuildIndices<N>{});
//	    }
//};


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
		//	...
		//	std::vector<data_t> samples;
		//	for(int i = 0; i < N; i++){
		//		val = priors[i].sample();
		//	}
		//}
};


// Variadic templates:
// template<typename... Arguments> class VariadicTemplate; or
// template<typename T, typename... Arguments> class VariadicTemplate; or
// template<typename... Arguments> void SampleFunction(Arguments..., ...)


// Vecchia idea:
// class DistributionWrapper{
// 	private:
// 		arma::vec parameters;
// 		shared_ptr<DistributionWrapper> prior = nullptr;
// 		arma::vec density() virtual;
// 	public:
// 		arma::vec sample_from_prior();
// 		Distribution(arma::vec parameters)
// }

template<state_t, data_t>
class DistributionWrapper{

	real eval(data_t x){
		...
	}
	void update(std::vector<state_t> vals){

	}
}


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


// Normal likelihoood, Normal Inverse Gamma Hierarchy
template<typename Hypers> //Hypers = TupleWrapper, distro, ...
class NNIGHierarchy {
protected:
	std::tuple<real, real> state;
	//std::tuple<> params; // mu_0, k, alpha, beta
	qualcosa che mi rappresenta i parametri ma possono cambiare
	Hypers* hypers;
public:
	double log_like(real datum) {
		return stan::math::normal_lpdf(datum, state-1, state-2);
	}

	void sample(std::vector<real> data) {
		// prendo mu_0, k0, alpha0, beta0 da hypers
		mu_0 = hypers. // current value
		// calcolo mu_post, k_post, alpha_post, beta_post 

		real sigmaNew = stan::math::inv_gamma_rng(alpha_post, beta_post, rng);
		real muNew = stan::math::normal_rng(mu_post, kpost * sigmaNew, rng);
		state(0) = muNew;
		state(1) = sigmaNew;
	}

	std::tuple getState();
}


template<Hierarchy, Hypers>
class Neal8 {
	int m=3;
	std::vector<Hierarchy> hierarchies;
	std::vector<real> data;
	std::array<n,int> allocations; // the c vector
	std::vector<data_t> unique_values;
	std::array<m,data_t> aux_unique_values;

	Hypers hypers;

	void sample_c {
		for (i=0; i<numData; i++) {
			old = c[i]
			probas = Eigen::VectorXd(numClusters);
			for (c = 0; c < numClusters; c++) {

				// meglio in logscale
				probas(c) = dataPerClusters[c] * hierarchies[c].like(data[i]);
			}

			probas;

			c[i] = stan::math::categorical_rng(probas, rng);
			dataPerClusters[old] -= 1;
			dataPerClusters[c[i]] += 1;
		}
	}

	void sample_clustervals(){
		for (c = 0)
	}

	void run(){
		initalize();
		while(!converged){
			step();
		}
	}

	void initalize(){
		//...
		return;
	}

	void step(){
		sample_allocations();
		sample_unique_values();
	}

	int unique_sans_i(int i){ // counts unique values by scrolling allocations
		// without the i-th unit
		// ...
		return k;
	}

	void sample_allocations(){
		for(int i=0; i<n; i++){ // for each data unit
			//data[i]
            std::vector<unsigned int> temp;
            // ideas: our own for loop for k and bool (ci is a singleton)
            // ideas: function from std count distinct values in vector
            // ideas: using a (multi)map?
            for(int j=0; j<n; j++){
                if(j == i)
                    continue;

            }

			int k = unique_sans_i(i);
			int h = k+m;
		}
	}
}

int main() {
	Neal8<NNIGHierarchy<HypersFixed>> sampler;
}


// Fai caso base: neal8 w/ normal-normalinvgamma, poi generalizza


// Modello N-NIG: kernel gaussiano, G0 N-IG:
// f ~ N(mu,sig^2)
// (mu,sig^2) ~ G
// G ~ DP(M, G0)  with G0 ~ N-IG


class Hypers{
    real current_value;
    void sample(){...} // prende a caso
}

