make_adjacency_full <- function(grid_size){
  grd <- .make_grid(grid_size, "full")
  .make_adjacency(grd)
}

make_adjacency_interior <- function(grid_size){
  grd <- .make_grid(grid_size, subset="int")
  .make_adjacency(grd)
}

make_adjacency_exterior <- function(grid_size){
  grd <- .make_grid(grid_size, subset="ext")
  .make_adjacency(grd)
}

.make_grid <- function(grid_size, subset=c('full', 'int','ext'))
{
  sfc = sf::st_sfc(sf::st_point(c(0,0)), sf::st_point(c(grid_size,grid_size)))
  grd =  sf::st_make_grid(sfc, n = c(grid_size,grid_size))
  
  
  if (subset != "full") {
    
    border <- get_border(grid_size)
    mask <- if(subset == "int") !border else border
    
    return(grd[mask])
  }
  
  grd
  
}

.make_adjacency <- function(x, .fun=st_rook, style="B"){
  
  sfc_nb = .fun(x)
  
  sf_nb = as_nb_sgbp(sfc_nb)
  W = spdep::nb2mat(sf_nb, style=style)
  
  W
}

get_border <- function(grid_size)
{
  border <- which((1:grid_size^2 %% grid_size == 1) | (1:grid_size^2 %% grid_size == 0))
  border <- c(border, 1:grid_size, (grid_size^2-grid_size+1):grid_size^2)
  
  
  1:grid_size^2 %in% border
  
}

st_rook = function(a, b = a) sf::st_relate(a, b, pattern = "F***1****")

as_nb_sgbp <- function(x, ...) {
  attrs <- attributes(x)
  x <- lapply(x, function(i) { if(length(i) == 0L) 0L else i } )
  attributes(x) <- attrs
  class(x) <- "nb"
  x
}