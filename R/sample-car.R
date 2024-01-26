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