class PrecMat
{
protected:
    LLT<Eigen::MatrixXd> cho_factor;
    Eigen::MatrixXd cho_factor_eval;
    double log_det;

public:
    PrecMat::PrecMat(const Eigen::MatrixXd &prec) : prec(prec) {
        cho_factor = LLT<Eigen::MatrixXd>(prec);
        cho_factor_eval = cho_factor.matrixL();
        const Eigen::VectorXd &diag = cho_factor_eval.diagonal();
        log_det = 2 * log(diag.array()).sum();
    }

    LLT<Eigen::MatrixXd> get_cho_factor() const {return cho_factor;}
    const Eigen::MatrixXd &get_cho_factor_eval() const {return cho_factor_eval;}
    double get_log_det() const {return log_det;}


    double multi_normal_lpdf(const std::vector<Eigen::VectorXd> &x,
        const Eigen::VectorXd &mu, const PrecMat &sigma) {
        int n = x.size();
        double out = sigma.get_log_det() * n;
    
        const Eigen::MatrixXd &cho_sigma = sigma.get_cho_factor_eval();
    
        std::vector<double> loglikes(n);
        for (int i = 0; i < n; i++)
        {
            loglikes[i] = (cho_sigma * (x[i] - mu)).squaredNorm();
        }
    
        out -= std::accumulate(loglikes.begin(), loglikes.end(), 0.0);
    
        return 0.5 * out;
    }

};
