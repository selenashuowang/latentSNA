#' Gibbs sampling of Theta
#'
#' A Gibbs sampler for updating the Person latent effect Theta.
#'
#' @usage rTheta(H, beta, Alpha, Theta ,U, Stheta, Sutheta, Su, s1)
#' @param H N X P normal behavior matrix
#' @param beta P X 1 behavior intercept vector
#' @param Alpha 1 X 1 intercept vector
#' @param Theta current value of Theta
#' @param U matrix containing current value of U
#' @param s1 behavior variance
#' @param Stheta covariance of Theta
#' @param Sutheta covariance between U and Theta
#' @param Su matrix containing covariance of U
#' @return \item{Theta}{a new value of Theta}
#' @author Selena Wang
#' @export rTheta_b
rTheta_b<-function(H, beta, Alpha, Theta ,U, Stheta, Sutheta, Su, s1)
{
  D<-nrow(Theta) ; n<-ncol(Theta); m = ncol(H)

  t<-H - matrix(rep(1,n)) %*% t(beta)


  ivTheta <- solve(Stheta - Sutheta %*% solve(Su) %*% t(Sutheta))


  invUTheta <- (-1)*solve(Stheta - Sutheta %*% solve(Su) %*% t(Sutheta)) %*% Sutheta %*% solve(Su)

  invThetaU <- (-1)* solve( Su - t(Sutheta) %*% solve(Stheta) %*% Sutheta) %*% t(Sutheta) %*% solve(Stheta)

if(D>1){
  ## update each Theta[i,]
  for(i in rep(sample(1:n),4))
  {
    l<-(rowSums(Alpha*t[i,], na.rm = TRUE) /s1 - .5*invUTheta %*% matrix(U[i,])
        - .5*t(invThetaU) %*% matrix(U[i,]))

    iQ<- solve(  ivTheta + ( Alpha %*% t(Alpha) )/s1  )

    Theta[,i]<- iQ%*%l + t(chol(iQ))%*%rnorm(D)

  }
}else{
  for(i in rep(sample(1:n),4))
  {
    l<-(sum(Alpha*t[i,], na.rm = TRUE) /s1 - .5*invUTheta %*% matrix(U[i,])
        - .5*t(invThetaU) %*% matrix(U[i,]))

    iQ<- solve(  ivTheta + sum( Alpha^2 )/s1  )

    Theta[,i]<- iQ%*%l + t(chol(iQ))%*%rnorm(D)

  }
}






  t(Theta)
}







