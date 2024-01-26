# Simulate a count matrix for analysis in spatial 
# Goals here are to load in "true parameters", generate nxn adjacency matrix
# Output a voyager object

sample_car_scaled <- function(W, rho, precision)
{
  dim_ <- dim(W)[1]
  
  # Rescale matricies
  D <-  diag(1/MatrixGenerics::rowSums(W))
  I <- diag(dim_)
  W <- W/MatrixGenerics::rowSums(W)
  
  .sample_mvnorm_prec(rep(0, dim_), precision*D %*% (I-rho*W))
}

sample_car <- function(W, rho, precision)
{
  dim_ <- dim(W)[1]
  D <- diag(MatrixGenerics::rowSums(W))
  
  .sample_mvnorm_prec(rep(0, dim_), precision * (D-rho*W))
}

.sample_mvnorm_prec <- function(mu, Omega)
{
  StanHeaders::stanFunction("multi_normal_prec_rng", mu=mu, S=Omega)
}

inv_logit <- function(x) 1/(1+exp(-x))

simulate_data <- function(
    grid_size, 
    n_celltypes, n_genes,
    rho, precision,
    empty_border = TRUE)
{
  # generate adjacency matrix and probability of capture
  message("Generate adjacency matrix...")
  W <- make_adjacency_full(grid_size)

  logit_prob_capture <- sample_car(W, rho, precision)
  message("Simulate technical parameter...")
  prob_capture <- inv_logit(logit_prob_capture)
  
  # meta  <- make_metadata(grid_size, n_celltypes, empty_border)
  message("Generate count matrix...")
  
  make_counts(
    grid_size=grid_size, 
    n_celltypes=n_celltypes, 
    n_genes=n_genes,
    prob_capture=prob_capture,
    empty_border=empty_border)
  
}
  
get_border <- function(grid_size)
{
  border <- which((1:grid_size^2 %% grid_size == 1) | (1:grid_size^2 %% grid_size == 0))
  border <- c(border, 1:grid_size, (grid_size^2-grid_size+1):grid_size^2)
 

  1:grid_size^2 %in% border
  
}

