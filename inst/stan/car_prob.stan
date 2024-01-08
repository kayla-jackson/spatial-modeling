data {
  int<lower=0> count_n; // number of fitted points
  array[count_n] int<lower=0> counts; // the data
  
  vector<lower=0, upper=1>[count_n] rates; // probability of capture
  
  real bs_mu;
  real mu_mu;
  real bs_var;
  real mu_var;
  
}

// The parameters accepted by the model
parameters {
  real<lower=0> mu; // mean to be estimated per gene
  real log_bs; // burst size to be estimated per gene
}

transformed parameters{
  real<lower=0> burst_size;
  burst_size = pow(10, log_bs);
}

model {
  counts ~ neg_binomial_2(mu, burst_size*rates);
  
  mu ~ normal(mu_mu, mu_var);
  log_bs ~ normal(bs_mu, bs_var);

}
