class Distribution{
	arma::vec parameters;
	double eval() virtual;
}

class GaussianDistribution: public Distribution{
    double eval() override{
			return 1/(2*pi)^0.5 * ...;
		}
}


//POI CREO OGGETTI GAUSSIAN, etc:
Distribution gaussian;


class DistributionWithHypers{
	std::vector<shared_ptr<DistributionWithFixed>> hyperparameters;
	arma::vec sample();
}

class DistributionWithFixed: public DistributionWithHypers{
	arma::vec parameters;
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
