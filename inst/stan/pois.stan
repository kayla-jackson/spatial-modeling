// use for case study NYC data subset
data {
  int<lower=0> N;
  int<lower=0> N_edges;
  // array[N_edges] int<lower=1, upper=N> node1; // node1[i] adjacent to node2[i]
  // array[N_edges] int<lower=1, upper=N> node2; // and node1[i] < node2[i]
  // real<lower=0> mu;
  vector[N] mu;
  matrix[N, N] D; 
  matrix[N, N] W;
  array[N] int<lower=0> counts; // count outcomes
}

transformed data{
  vector[N] zeros; 
  for (i in 1:N){
    zeros[i] = 0;
  }
}
parameters {
  //real beta0; // intercept
  vector[N] phi; // spatial effects
  real<lower=0> car_scale;
  real<lower=0,upper=1> car_rho;
}

transformed parameters{
  vector[N] rates; 
  rates = 1 / (1+exp(-phi));
}

model {
  counts ~ poisson(mu .* rates);
  //beta0 ~ normal(0, 1);
  
  //prior for phi
  //target += -0.5 * dot_self(phi[node1] - phi[node2]);
  car_scale ~ normal(0,1);
  car_rho ~ beta(3, 1);
  
  phi ~ multi_normal_prec(zeros, car_scale * (D-car_rho*W));
  
  // soft sum-to-zero constraint on phi)
  // sum(phi) ~ normal(0, 0.001 * N); // equivalent to mean(phi) ~ normal(0,0.001)
}

generated quantities {
  
  //vector[N] eta = phi;
  vector[N] prob_ppc = 1/(1+exp(-phi));
  array[N] int<lower=0> counts_ppc = poisson_rng(mu .* prob_ppc);
  
}

