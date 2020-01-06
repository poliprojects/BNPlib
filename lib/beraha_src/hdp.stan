functions {
  real gammaMix_log(real x, real w, real a1, real b1, real a2, real b2) {
    return log_mix(w,gamma_lpdf(x | a1, b1),gamma_lpdf(x | a2, b2));
  }
}

data {
  int<lower=0> H; // number of components in the DP
  int<lower=0> J; // number of groups
  int<lower=0> max_num_samples;
  int<lower=0> num_samples[J]; // number of samples in each group
  matrix[J, max_num_samples] samples;
  real a1;
  real b1;
  real a2;
  real b2;
}

parameters {
  real<lower=0, upper=1> w;
  real<lower=0.2> alpha_0;
  vector<lower=0, upper=1>[H-1] nus[J];
  vector<lower=0, upper=1>[H-1] nu_top;
  vector[H] means; // cluster means
}

transformed parameters {
  simplex[H] weights[J];
  simplex[H] weights_top;
  real prod1_nu = 0;

  for (j in 1:J) {
    weights[j][1] = nus[j][1];
    prod1_nu = 1 - nus[j][1];
    for (h in 2:(H-1)) {
      weights[j][h] = nus[j][h] * prod1_nu;
      prod1_nu *= (1 - nus[j][h]);
    }
    weights[j][H] = fmax(0.0, 1 - sum(weights[j][1:(H-1)]));
  }

  weights_top[1] = nu_top[1];
  prod1_nu = 1 - nu_top[1];
  for (h in 2:(H-1)) {
    weights_top[h] = nu_top[h] * prod1_nu;
    prod1_nu *= (1 - nu_top[h]);
  }
  weights_top[H] = fmax(0.0, 1 - sum(weights_top[1:(H-1)]));
}

model {

  // hyperparams for alpha_0 ~ w*Gamma(a1, b1) + w*Gamma(a2, b2)

  // hyperparams for G0 | gamma, H
  // H ~ N(0, 5);
  real sigmaH = 1;
  real gamma = 2.0;
  // real alpha_0 = 1.0;

  // Top level DP
  for (h in 1:H) {
    means[h] ~ normal(3, sigmaH);
  }

  nu_top ~ beta(1, gamma);

  w ~ uniform(0, 1);
  alpha_0 ~ gammaMix(w, a1, b1, a2, b2);

  // Bottom level DPs
  for (j in 1:J) {
    nus[j] ~ beta(1, alpha_0);
  }

  for (j in 1:J) {
    for (i in 1:num_samples[j]) {
      real partial_sums[H];
      for (h in 1:H) {
        partial_sums[h] = log(weights[j][h]) + normal_lpdf(samples[j, i] | means[h], 1.0);
      }
      target += log_sum_exp(partial_sums);
    }
  }
}
