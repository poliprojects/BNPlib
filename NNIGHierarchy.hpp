// Normal likelihoood, Normal Inverse Gamma hierarchy
template<class Hypers> //Hypers = TupleWrapper, distro, ...
class NNIGHierarchy {
protected:
    using state_t= std::array<par_t,2>;

    std::mt19937 rng;

    state_t state; // current values for F's parameters: mu, sigma

    std::shared_ptr<Hypers> hypers; // current values for G0's parameters:
                                    // mu_0,Lambda0, alpha, beta

public:
    // Constructors:
    ~NNIGHierarchy() = default;
    NNIGHierarchy(std::shared_ptr<Hypers> hypers): hypers(hypers) {}

    // Getters/setters:
    state_t get_state(){return state;}
    void set_state(const state_t &s){state = s;}
    void set_state(int pos, par_t val){state[pos] = val;}

    int get_count(){return hypers.use_count();}


    double log_like(data_t datum) {
        return exp(stan::math::normal_lpdf(datum, state[0], state[1]));
    }

    void draw() {
        float sigmaNew = stan::math::inv_gamma_rng(hypers->get_alpha0(),
            hypers->get_beta0(), rng);
        float muNew = stan::math::normal_rng(hypers->get_mu0(),
            sigmaNew/hypers->get_lambda(), rng);
        state[0] = muNew;
        state[1] = sigmaNew;
    }


    void sample_given_data(std::vector<data_t> data) {
        // Get current values of parameters
        auto mu0    = hypers->get_mu0();
        auto Lambda0  = hypers->get_lambda();
        auto alpha0 = hypers->get_alpha0();
        auto beta0  = hypers->get_beta0();

        arma::vec temp = normalGammaUpdate(data, mu0, alpha0, beta0, Lambda0);

        auto mu_post = temp(0);
        auto alpha_post = temp(1);
        auto beta_post = temp(2);
        auto postLambda = temp(3);

        // Get a sample
        par_t sigma_new = stan::math::inv_gamma_rng(alpha_post, beta_post, rng);
        par_t mu_new = stan::math::normal_rng(mu_post, sigma_new/postLambda,
            rng); //? is it ok /postLambda?
        state[0] = mu_new;
        state[1] = sigma_new;
    }



    arma::vec normalGammaUpdate(arma::vec data, double mu0, double alpha0,
        double beta0,double Lambda0) {
            double mu_post, alpha_post, beta_post, postLambda;
            int n = data.size();
            if (n == 0)
                return arma::vec{mu0, alpha0, beta0, Lambda0};
    double ybar = arma::mean(data);
    mu_post = (Lambda0 * mu0 + n * ybar) / (Lambda0 + n);
    alpha_post = 1.0 * alpha0 + 1.0 * n / 2;

    // arma::var(x, 1) divides by n, not n-1
    double ss = n * arma::var(data, 1);

    beta_post = (beta0 + 0.5 * ss + 0.5 * Lambda0 /
        (n + Lambda0) * n * std::pow((ybar - mu0), 2));

    postLambda = Lambda0 + n;

    return arma::vec{mu_post, alpha_post, beta_post, postLambda};
  }

};
