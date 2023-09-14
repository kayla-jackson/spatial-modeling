# Pull more than one variable from a data frame

pull2 <- function(data, vars = -1){
  
  data <- fortify_data(data)
  vars <- get_vars_exist(names = names(data), !!enquo(vars))
  
  # only length of existing vars
  n <- length(vars)
  if (n == 0L) {
    stop("no variables found")
  }
  
  if (n == 1L){
    return(data[[vars]])
  }
  
 purrr::map(vars, ~ data[[.x]])
  
  
}

get_vars_exist <- function(names, vars = -1L,..., env = current_env()){
  
  check_dots_empty()
  
  expr_ <- enquo(vars)
  if (quo_is_missing(expr_)) {
    stop("must supply variables")
  }
  
  loc <- eval_tidy(expr_, set_names(seq_along(names), names))
  # print(loc)
  if (any(loc < 0L)) {
    n <- length(names)
    where <- which(loc < 0L)
    
    sub_loc <- loc[where]
    
    loc[where] <- n + sub_loc + 1L
  }
  
  names[loc] 
}


fortify_data <- function(x) {
  data.frame(x)
}
