#' Simulation from a Wishart distribution
#' 
#' Simulates a random Wishart-distributed matrix
#' 
#' 
#' @usage rwish(S0, nu = dim(S0)[1] + 2)
#' @param S0 a positive definite matrix
#' @param nu a positive integer
#' @return a positive definite matrix
#' @author Selena Wang
#' 
#' 
#' 
#' 
rwish <-
function(S0,nu=dim(S0)[1]+2)
{
  # sample from a Wishart distribution 
  # with expected value nu*S0 
  sS0<-chol(S0)
  Z<-matrix(rnorm(nu*dim(S0)[1]),nu,dim(S0)[1])%*%sS0
  t(Z)%*%Z
}
