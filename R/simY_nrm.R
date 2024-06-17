#' Simulate a normal behavior matrix
#' 
#' Simulates a normal behavior matrix
#' 
#' 
#' @usage simY_nrm(EY, s1)
#' @param EY matrix giving the expected value of the behavior matrix
#' @param s1 variance
#' @return a N by P matrix
#' @author Selena Wang
#' @export simY_nrm
simY_nrm <-
function(EY,s1) 
{
  YS<-simH(EY,s1) 
  YS
}
