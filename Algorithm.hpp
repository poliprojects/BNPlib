// GENERIC ALGORITHM
template<MixingMode? M, Hierarchy H>
class Algorithm{
	private:
		unsigned int maxiter;
		void step(unsigned int iter) {
			std::cout << "Step " << iter << std::endl;
			sample_clusters(); // c_i, allocation
			sample_parameters(); // phi_c
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
