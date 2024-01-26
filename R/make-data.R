.params <- readr::read_csv("./data/parameters.csv")

#' @param grid_size Used to create a regular lattice of dimension grid_size x grid_size. 
#'        If `empty_border = FALSE`, then the total number of 'cells' in the
#'        simulated data will be equal to the square of `grid_size`
#' @param n_celltypes Number of cell types to simulate
#' @param n_genes Number of genes to simulate
#' @param empty_border Does the border of the lattice represent 'empty' spots?
make_metadata <- function(grid_size, n_celltypes, empty_border=TRUE)
{
  
  border <- NULL
  n_cells <- grid_size^2
  
  if (empty_border) {
    border <- get_border(grid_size)
    n_cells <- grid_size^2 - sum(border)
    
  }
  
  cellids <- label_cellids(grid_size)
  types <- label_celltypes(n_celltypes, n_cells, border_i=border)
  
  data.frame(cellid = cellids, type = types)
}

make_genedata <- function(
    grid_size, 
    n_celltypes, 
    n_genes, 
    prob_capture, 
    empty_border=TRUE)
{
  # eventually, call: make_data_("gene", ...)
  # add checks on req'd args...
  
  if (empty_border) border = get_border(grid_size) else border = NULL
 
  
  ct_labs <- label_celltypes(grid_size, n_celltypes, border_i=border)
  ct_obs <- unique(ct_labs[which(ct_labs != "border")])
  n_celltypes_obs <- length(ct_obs)
  
  dimnames <- list(label_genes(n_genes), label_cellids(grid_size))
  
  cell_means <- matrix(0, nrow=n_genes, ncol=grid_size^2, dimnames=dimnames)
  cell_disp  <- matrix(1, nrow=n_genes, ncol=grid_size^2, dimnames=dimnames)
  
  resample <- nrow(.params) < n_genes
  
  for (type in ct_obs) {
    
    mask <- which(ct_labs == type)
    rows <- sample(nrow(.params), n_genes, replace=resample)
  
    means <- .params[[2]][rows]
    bs    <- .params[[3]][rows] 
    
    cell_means[,mask] <- means
    cell_disp[,mask]  <- bs
  
  }
  # Update bs to include probability of capture.
  # If border cell, bs = prob_capture
  # Not most efficient. Fix later
  cell_disp <- t(t(cell_disp) * prob_capture)
  
  # Estimate cell means for the border, simply the average of the neighbors
  if (!is.null(border)){
    W <- make_adjacency_full(grid_size)
    means <- crossprod(W, t(cell_means))/MatrixGenerics::rowSums(W)
    
    cell_means[,border] <- t(means[border,])
  }
  
  
  list(means=cell_means, disp=cell_disp)
  
}

make_counts <- function(
    grid_size, 
    n_celltypes, 
    n_genes, 
    prob_capture, 
    empty_border=TRUE)
{
  
  gene.params <- make_genedata(grid_size, n_celltypes, n_genes, prob_capture, empty_border)
  
  nb.counts <- matrix(
    rnbinom(grid_size^2*n_genes, size=gene.params$means, prob=1/(1+gene.params$disp)),
    nrow=n_genes, ncol=grid_size^2)
  
  if(empty_border){
    # Sample counts along the border
    border <- get_border(grid_size)
    pois.counts <- matrix(
      rpois(grid_size^2*n_genes, lambda=gene.params$means*gene.params$disp),
      nrow=n_genes, ncol=grid_size^2
    )

  }
  
  nb.counts[,border] <- pois.counts[,border]
  dimnames(nb.counts) <- dimnames(gene.params$means)
  
  message("Returning counts...")
  nb.counts
}

