


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

transformed data{
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
  real<lower=0> mu; // mean to be estimated per gene
  real log_bs; // burst size to be estimated per gene
  vector[N] phi;
}

transformed parameters{
  real<lower=0> burst_size;
  vector[count_n] rates;
  vector[count_n] phi_sub;
  
  burst_size = pow(10, log_bs);
  phi_sub = phi[count_inds];
  
  rates = 1 ./ (1+exp(-phi_sub));
}

model {
  counts ~ neg_binomial_2(mu, burst_size*rates);
  
  mu ~ normal(mu_mu, mu_var);
  log_bs ~ normal(bs_mu, bs_var);
  
  phi ~ multi_normal_prec(zeros, car_scale * (D-car_rho*W));
}
