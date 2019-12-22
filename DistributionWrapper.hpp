namespace bnplib{
	using data_t = double;
	using param_t = std::vector<double>;
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

	//qualcosa che mi rappresenta i parametri ma possono cambiare
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
	std::vector<param_t> unique_values;
	std::array<m,param_t> aux_unique_values;

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

		// initialize cardinalities
		std::vector<unsigned int> card(unique_values.size(), 0); // cardinalities of unique values

		for(int j=0; j<n; j++){
			card[allocations[j]]+=1;
		}


		int k, start, n_unique;

		for(int i=0; i<n; i++){ // for each data unit

			//data[i]
			// ideas: our own for loop for k and bool (ci is a singleton)
			// ideas: function from std count distinct values in vector
			// ideas: using a (multi)map?
			start=0;
			n_unique=unique_values.size();
			if(card[allocations[i]]==1){ // datum i is a singleton
			k=n_unique-1;
			aux_unique_values[0]=unique_values[allocations[i]];
			singleton=1;}
			else k=n_unique;
			card[allocations[i]]-=1;

			for(int j=singleton;j<m; j++ ){
			aux_unique_values[j][0] = stan::math::normal_rng(hypers.get_m0(), hypers.get_k0(), rng);
			aux_unique_values[j][1] = stan::math::inv_gamma_rng(hypers.get_alpha0(),hypers.get_beta0(), rng);
			}
			probas = Eigen::VectorXd(k+m); //			std::vector<double> probas(k+m,0);

			for(int k=0; k<n_unique ; k++){ // if datum i is a singleton card[k] when k=allocations[i] is equal to 0 ->probas[k]=0
				probas(k)=card[k]*stan::math::normal_cdf(data[i],unique_values[k][0],unique_values[k][1])/(n-1+ mixture.get_totalmass());
			}
			for(int k=0; k<m ; k++){
				probas(n_unique+k)=(mixture.get_totalmass()/m)*stan::math::normal_pdf(data[i],unique_values[k][0],unique_values[k][1])/(n-1+ mixture.get_totalmass());
			}

			int c_new=  stan::math::categorical_rng(probas,rng) //std::discrete_distribution(probas); // in stats??
			if(singleton){

				if (c_new>=n_unique){ // SINGLETON - AUX
					unique_values[allocations[i]]=aux_unique_values[c_new-n_unique];
					card[allocations[i]]+=1;
				}
				else // SINGLETON - OLD VALUE
		    {unique_values.erase(unique_values.begin()+allocations[i]-1);
				card.erase(card.begin()+allocations[i]-1);
				card[c_new]+=1;

				int tmp=allocations[i];
				allocations[i]=c_new;
				for(auto &c:allocations){ // relabeling
					if (c>tmp)
					c-=1;
				}
			}

			}
			else{

				if (c_new>=n_unique){ // NOT SINGLETON - AUX
					unique_values.push_back(aux_unique_values[c_new-n_unique]);
					card.push_back(1);
					allocations[i]=n_unique+1;
			}
			else{ // NOT SINGLETON - OLD VALUES
				allocations[i]=c_new;
				card[c_new]+=1;
			}

		}
			}
		} // end sample allocation


void sample_unique_values(){
std::vector<std::vector<unsigned int>> clust_idxs;
for(int i=0; i<n; i++) // save different cluster in each row
	clust_idxs[allocations[i]].push_back(i);


for (auto &row: clust_idxs) {
	std::vector<data_t> curr_data;
	for (auto &idx: row){
		curr_data.push_back(data[idx]);
	}
	// draw PHI_c given curr_data;
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
