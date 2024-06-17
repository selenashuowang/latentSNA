#' Gibbs sampling of behaivor parameters Xi
#'
#' A Gibbs sampler for updating the behavior parameters.
#'
#' @usage rXi(H, beta, Alpha, Theta, mxi, Sigmaxi, s1 = 1)
#' @param H normal behavior matrix
#' @param beta  behavior intercept vector
#' @param Alpha  behavior intercept vector
#' @param Theta  current value of Theta
#' @param mxi  vector of prior for the mean of Xi
#' @param Sigmaxi matrix of prior for the variance of Xi
#' @param s1 behavior variance
#' @return \item{beta}{a new value of beta} \item{Alpha}{a new value of Alpha}
#' @author Selena Wang
#' @export rXi
rXi<-function(H, beta, Alpha, Theta, mxi, Sigmaxi, s1)
{
  D<-ncol(Theta) ; n<-nrow(Theta); m = ncol(H)

  ## G is a N X (D+1) matrix
  G= cbind(Theta,rep(1,n))
  ## Xi is a (D+1) X M matrix

  H[is.na(H)]=0

  Xi= rbind(t(Alpha), t(beta))
  ## update each Xi[i,]
  for(i in rep(sample(1:m),4))
  {
    l<- t(G) %*% matrix(H[,i], nrow = n)/s1 +  solve(Sigmaxi) %*% mxi
    iQ<- solve( ( solve(Sigmaxi) + ( t(G) %*% G )/s1 ) )

    # a=G[1,] %*% t(G[1,])
    # for(i in 2:N){
    #  a=a+G[i,] %*% t(G[i,])
    # }
    # tmp=iQ%*%l + t(chol(iQ))%*%rnorm(D+1)
    # A=c(tmp[1:D,1]>0,TRUE)
    # Xi[,i][A] <- tmp[A]

    # tmp=iQ%*%l + t(chol(iQ))%*%rnorm(D+1)
    # Xi[,i][1:D][tmp[1:D]>0]=tmp[1:D][tmp[1:D]>0]
    #
    # Xi[,i][D+1] = tmp[D+1]

    Xi[,i] = iQ%*%l + t(chol(iQ))%*%rnorm(D+1)

  }




  if(D==1){
    return(list(beta=matrix(Xi[nrow(Xi),]), Alpha=matrix(Xi[1:D,])))
  }
  if(D>1){
    return(list(beta=matrix(Xi[nrow(Xi),]), Alpha=t(Xi[1:D,])))
  }



}







