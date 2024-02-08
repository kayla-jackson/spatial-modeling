data {
  int<lower=0> N; // number of areal units
  int<lower=0> count_n; // number of fitted points
  array[count_n] int<lower=0> counts; // the data
  array[count_n] int<lower=0, upper=N> count_inds; //location of the counts?

  
  real bs_mu;
  real mu_mu;
  real bs_var;
  real mu_var;
  
  real car_scale;
  real car_rho;
  matrix[N, N] W;
}

transformed data {
  vector[N] zeros;
  matrix<lower = 0>[N, N] D;
  {
    vector[N] W_rowsums;
    for (i in 1:N) {
      W_rowsums[i] = sum(W[i, ]);
    }
    D = diag_matrix(W_rowsums);
  }
  zeros = rep_vector(0, N);
}


// The parameters accepted by the model
parameters {
  real log_mu; // relative transcription rate per gene
  real log_bs; // burst size to be estimated per gene
  vector[N] phi;
}

transformed parameters {
  real<lower=0> beta_ = pow(10, log_bs);
  real<lower=0> trans_rate = pow(10, log_mu);
  
  vector[count_n] phi_sub = phi[count_inds];
  
  vector[count_n] rates;
  rates = 1 / (1+exp(-phi_sub));
}

model {
  for (n in 1:count_n)
    counts[n] ~ neg_binomial(trans_rate, beta_*rates[n]);

  log_mu ~ normal(mu_mu, mu_var);
  log_bs ~ normal(bs_mu, bs_var);
  
  phi ~ multi_normal_prec(zeros, car_scale * (D-car_rho*W));
}
