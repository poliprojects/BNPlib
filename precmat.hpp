class PrecMat
{
protected:
    LLT<MatrixXd> cho_factor;
    MatrixXd cho_factor_eval;
    double log_det;

public:
    PrecMat::PrecMat(const MatrixXd &prec) : prec(prec) {
	    cho_factor = LLT<MatrixXd>(prec);
   		cho_factor_eval = cho_factor.matrixL();
    	const VectorXd &diag = cho_factor_eval.diagonal();
    	log_det = 2 * log(diag.array()).sum();
	}

    LLT<MatrixXd> get_cho_factor() const {return cho_factor;}
    const MatrixXd &get_cho_factor_eval() const {return cho_factor_eval;}
    double get_log_det() const {return log_det;}
};
