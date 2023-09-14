# Some helper functions
matcher <- function(x) names(ct_sf)[-1][as.integer(x)]

is_zero <- function(x)  x == 0

mse <- function(x, y=0){
  mean((x-y)^2)
}

semivar <- function(gene, .inds, obj) {
  
  err <- map_dbl(.inds, \(x) {
    
    mse(
      counts(obj)[gene,][ x[[1]] ], 
      counts(obj)[gene,][ x[[2]] ]
    )
  })
  
  err/2
}

semivar2 <- function(feature, .inds, obj) {
  
  err <- map_dbl(.inds, \(x) {
    
    mse(
      colData(obj)[[feature]][ x[[1]] ], 
      colData(obj)[[feature]][ x[[2]] ]
    )
  })
  
  err/2
}
