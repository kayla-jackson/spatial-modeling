grid_size <- 10
n_celltypes <- 2
n_genes <- 100
rho <- 0
precision <- 1

pois.file <- fs::path("./inst/stan/pois", ext="stan")
car.file  <- fs::path("./inst/stan/car",ext="stan")

# Compile stan models
car.mod <- cmdstanr::cmdstan_model(car.file)
pois.mod <- cmdstanr::cmdstan_model(pois.file)

counts <- simulate_data(
    grid_size, 
    n_celltypes, n_genes,
    rho, precision,
    empty_border = TRUE
)

# Set up data for fitting poisson model
border <- get_border(10)
W <- make_adjacency_full(10)

inds <- rep(0,100)
for(i in seq(n_genes)){
  vals <- crossprod(W, counts[i,])/rowSums(W)
  if(all(vals[border,1] > 0)){inds[i] <- i}
}

inds <- which(as.logical(inds))
gene <- inds[2]
moran(counts[gene, border], col.w.ext, n=36, Szero(col.w.ext))

border.means <- crossprod(W, counts[gene,])/rowSums(W)
border.means <- border.means[border,1]

moran(border.means, col.w.ext, n=36, Szero(col.w.ext))

# Tests
grid_size <- 20
ps <- inv_logit(sample_car(make_adjacency_full(grid_size), 0.72, 1))
#cts <- rpois(36, ps[border]*border.means)
cts <- rpois(grid_size^2, ps*30)

pois.data <- list(
  N = grid_size^2, 
  count_n = grid_size^2, 
  count_inds = seq(grid_size^2),
  W = as.matrix(make_adjacency_full(grid_size)),
  mu = rep(mean(cts)/0.5, grid_size^2),
  counts = cts
)


# Fit Poisson model to empty spots

#fit.pois <- pois.mod$optimize(pois.data)
fit.pois <- pois.mod$sample(pois.data, chains=2, iter_sampling=1000)
fit.pois$summary(c("car_rho", 'car_rho_trans'))



mods <- lapply(inds, function(i){
  #border.means <- crossprod(W, counts[i,])/rowSums(W)
  #border.means <- border.means[border,1]
  pois.data$mu <-  rep(mean(counts[gene,]), sum(border))
  pois.data$counts <- counts[i, border]
  
  mod <- pois.mod$optimize(pois.data)
  mod
})


draws <- posterior::as_draws_matrix(fit.pois)
#yrep <- draws[, glue::glue('rates[{1:100}]')]
yrep <- draws[, c('car_rho', 'car_rho_trans')]
p <- data.frame(yrep) %>% 
  tibble::rownames_to_column('draw') %>%
  tidyr::pivot_longer(
    cols = contains("car_rho"), 
    names_to="param", 
    values_to="value") %>% 
  # mutate(
  #   ind = sub("prob_ppc.","", ind), 
  #   ind =  as.factor(sub("\\.","", ind))
  # ) %>%
  ggplot(aes(x = value)) + 
  geom_histogram(stat='density') +
  facet_wrap(~param)
  #ggplot(aes(ind, prob_val)) + 
  #geom_boxplot()

p +  geom_point(
  data=data.frame(ind = factor(seq(100), levels=seq(100)), prob_val = ps), 
  color = 'red')

y <- counts[10, border]
bayesplot::ppc_dens_overlay(y, yrep[1:200, ])

# Fit NB model to "cells"
nb.data <- 
  list(
    N = 100,
    W = W,
    count_n = length(!border),
    count_inds = which(!border),
    car_rho=hyperpars[1],
    car_scale=hyperpars[2],
    bs_mu = 0, 
    bs_var = 1, 
    mu_mu = 0, 
    mu_var = 5, 
    counts = y
  )
