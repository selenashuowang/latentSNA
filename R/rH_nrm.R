#' Simulate missing values in a normal behavior model
#' 
#' Simulates missing values under a behavior model
#' 
#' 
#' @usage rH_nrm(H, EH,s1, Y)
#' @param H a matrix, the current value of H
#' @param EH expected value of H
#' @param s1 behavior variance
#' @param Y behavior matrix
#' @return a behavior matrix, equal to  at non-missing values
#' @author Selena Wang
#' @export rH_nrm
rH_nrm<-function(H,EH,s1,Y)
{
  HS<-simY_nrm(EH,s1)
  H[is.na(Y)]<-HS[is.na(Y)]  # this isn't quite right if there is asymmetric missingness. 
  H
}

