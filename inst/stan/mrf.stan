// Trying to simulate MRF
functions {
  real icar_normal_lpdf(vector alphas, int N, int[] node1, int[] node2) {
    return -0.5 * dot_self(alphas[node1] - alphas[node2]);
 }
}

data {
  int N;
  int N_edges;
  int<lower=1, upper=N> node1[N_edges];
  int<lower=1, upper=N> node2[N_edges];
  
  matrix<lower = 0>[N, N] W;
  matrix<lower = 0>[N, N] D;
}

transformed data{
  vector[N] zeros;
  zeros = rep_vector(0, N);
}

parameters{
  vector<lower=0, upper=1>[N] theta;
  vector<lower=0>[N] alphas;
  vector<lower=0>[N] betas;
}

model{
  alphas ~ icar_normal_lpdf(N, node1, node2);
  betas ~  icar_normal_lpdf(N, node1, node2);
  
  //sum(alphas) ~ normal(0, 0.001);
  //sum(betas) ~ normal(0, 0.001);
  
  theta ~ beta(alphas, betas);
}
