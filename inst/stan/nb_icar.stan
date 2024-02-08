functions {
  real icar_normal_lpdf(vector phi, int N, array[] int node1,
                        array[] int node2) {
    return -0.5 * dot_self(phi[node1] - phi[node2]);
  }
}
data {
  int<lower=0> N; // number of areal units
  int<lower=0> N_edges;
  array[N_edges] int<lower=1, upper=N> node1; // node1[i] adjacent to node2[i]
  array[N_edges] int<lower=1, upper=N> node2; // and node1[i] < node2[i]
  
  array[N] int<lower=0> y; // count outcomes
  
  //hyperparameters
  real mu_var;
  real mu_mu;
  real bs_mu;
  real bs_var;
}

parameters {
  real<lower=0> mu; 
  real log_bs; // base 10 log of burst size
  real<lower=0> sigma;
  vector[N] phi; // spatial effects
}

transformed parameters {
  real<lower=0> burst_size;
  vector[N] rates;
  
  burst_size = pow(10, log_bs);
  
  vector[N] eta = sigma^2 * phi;
  rates = inv(1+exp(-eta));
}

model {
  for (n in 1:N)
    y ~  neg_binomial(mu, burst_size*rates[n]);
  
  mu ~ normal(mu_mu, mu_var);
  log_bs ~ normal(bs_mu, bs_var);
  
  sigma ~ normal(0.0, 3.0);
  phi ~ icar_normal(N, node1, node2);
  
  // soft sum-to-zero constraint on phi
  sum(phi) ~ normal(0, 0.001 * N);
}