label_celltypes <- function(grid_size, n_celltypes, border_i)
{
  # need checks on border_i? 
  if (!rlang::is_null(border_i)) {
    n_cells <- grid_size^2 - sum(border_i) 
  } 
  else {
    n_cells <- grid_size^2
  }
    
  if (n_celltypes > n_cells) n_celltypes <- n_cells
  
  rem <- n_cells %% n_celltypes
  ct_size <- rep(floor(n_cells/n_celltypes), n_celltypes)
  
  if(rem > 0L) ct_size[1:rem] <- ct_size[1:rem] + 1
  
  if (n_celltypes <= 26) {
    labels <- letters
    
  } 
  
  else {
    lab_len <- ceiling(n_celltypes/26)
    labels <- tidyr::expand_grid(t1=letters,t2=letters[1:lab_len])
    labels <- tidyr::unite(ct_labs, "lab", t1, t2, sep="")$lab
    
  }
  
  labels <- rep(labels[1:n_celltypes], ct_size)
  
  types <- rep("border", n_cells + sum(border_i))
  types <- if(length(border_i != 0L)) replace(types, !border_i, labels) else labels
  
  types
  
}

label_cellids <- function(grid_size){
  ids <- stringr::str_pad(1:grid_size^2, nchar(grid_size^2), pad="0")
  paste0("cell_", ids)
}

label_genes <- function(n_genes){
  ids <- stringr::str_pad(1:n_genes, nchar(n_genes), pad="0")
  paste0("gene_", ids)
}