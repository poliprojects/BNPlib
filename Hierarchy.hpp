class Distribution{
	arma::vec parameters;
	double eval();
}


class DistributionWithHypers: public Distribution{
	std::vector<shared_ptr<Distribution>> hyperparameters;
	arma::vec update();
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
