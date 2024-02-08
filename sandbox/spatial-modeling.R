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
  


