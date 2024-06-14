#' Gibbs update for latent effects covariance
#' 
#' @usage rSu(U,Su0=NULL,etau=NULL) 
#' @param U latent connectivity and behavior
#' @param Su0 prior (inverse) scale matrix for the prior distribution
#' @param etau prior degrees of freedom for the prior distribution
#' @author Selena Wang
#' @export rSu
#'
rSu<-function(U, Su0=NULL,etau=NULL) 
{
  
  Q=U
  if(is.null(Su0)){ Su0<-diag(ncol(Q))  } 
  if(is.null(etau)){ etau<-ncol(Q)+2 }

  S=solve(rwish(solve(etau*Su0+t(Q) %*% Q), etau+nrow(Q)))
  
    list("Su" = S)


}


