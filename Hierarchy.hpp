class MixingMode?{
	...
}

class Hierarchy{
	...
}

template<MixingMode? M, Hierarchy H>
class Algorithm{
	private:
		unsigned int maxiter;
		void step(unsigned int iter) {
			std::cout << "Step " << iter << std::endl;
			sample_clusters();
			sample_parameters();
			if(H->G0.has_prior())
				update_hyperparams();
			if(M.alpha_has_prior())
				update_alpha();
			if(iter > burnin)
				save();
		}
		void sample_clusters(){...}
		// etc...

	public:
		// Constructors...

		void run(){
			initialize();
			unsigned int iter = 1;
			while(!converged){
				step(iter);
				iter++;
				converged = ...;
			}
		}
}







class params {
 protected:
  double mu0, a, b, lambda; // state
  double mu0mean, mu0var, lambdaA, lambdaB;
  stats::rand_engine_t& engine = RandomEngine::Instance().get();
public:
  friend class NormalConjugateHierarchy;
  void update(arma::mat clusterVals);
  double getMu0() {return mu0;}
  double getA() {return a;}
  double getB() {return b;}
  double getLambda() {return lambda;}
};


class Hierarchy{
	DistributionWithHypers base_measure;
	DistributionWithHypers kernel;
}
