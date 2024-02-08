data {
  int<lower=0> N;
  int<lower=0> count_n; // number of fitted points
  array[count_n] int<lower=0, upper=N> count_inds; //location of the counts?

  vector[count_n] mu;
  matrix[N, N] W;
  array[count_n] int<lower=0> counts; // the data

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

parameters {
  vector[N] phi; // spatial effects
  real<lower=0> log10_car_scale;
  // real<lower=0> car_scale;
  real<lower=0,upper=1> car_rho;
}

transformed parameters{
  vector[count_n] rates;
  vector[count_n] phi_sub;
  
  real<lower=0> car_scale = pow(10, log10_car_scale);
  real<lower=-1> car_rho_trans = (2*car_rho) - 1;
  
  phi_sub = phi[count_inds];
  
  rates = 1 / (1+exp(-phi_sub));
  //rates = exp(phi);
  
}

model {
  
  car_rho ~ beta(1.6, 2.6);
  log10_car_scale ~ normal(0, 1);
  phi ~ multi_normal_prec(zeros, (D - car_rho_trans*W));
  
  for (n in 1:count_n)
    counts[n] ~ poisson(mu[n] * rates[n]);
  
}

generated quantities {

  vector[count_n] prob_ppc = 1/(1+exp(-phi_sub));
  array[count_n] int<lower=0> counts_ppc = poisson_rng(mu .* prob_ppc);

}

