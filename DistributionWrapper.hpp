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
// 		bool has_prior;
// 		shared_ptr<Distribution> prior = nullptr;
// 		arma::vec density() virtual;
// 	public:
// 		arma::vec sample_from_prior();
// 		Distribution(arma::vec parameters)
// }
