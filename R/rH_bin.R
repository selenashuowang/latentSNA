#' Simulate H 
#' 
#' Simulates a random latent matrix H given its expectation
#' and a behavior matrix Y
#' 
#' 
#' @usage rH_bin(H,EH,Y,s1=1)
#' @param H a n X m matrix, the current value of H
#' @param EH expected value of H
#' @param Y n X m binary item response matrix
#' @param s1 item response variance
#' @return a n X m matrix, the new value of H
#' @author Selena Wang
#' @export rH_bin
rH_bin <- function(H,EH,Y,s1=1 )
  { 

    N <- nrow(Y); M <- ncol(Y)
    E<-matrix(rnorm(N*M, sd = sqrt(s1)),N,M)
    HP<-EH + E
    

    A<-(( sign(HP) == sign(Y-.5)) & !is.na(Y)) 
    A[is.na(Y)]<-TRUE
    H[A]<-HP[A]
    

    H
  }

