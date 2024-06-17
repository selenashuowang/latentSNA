#' Gibbs update for behavior variance
#' 
#' Gibbs update for behavior variance
#' 
#' 
#' @usage rs1(Tm, offset=0,nu1=NULL,s10=NULL)
#' @param Tm a list of V X P normal behavior matrix
#' @param nu1 prior degrees of freedom 
#' @param s10 prior estimate of s1
#' @return a new value of s1
#' @author Selena Wang
#' @export rs1_b
rs1_b <-function(Tm,offset = offset, nu1=NULL,s10=NULL)
{ 
  
  tmp=Tm-offset
  tmp.1=sum(as.numeric(tmp)^2,na.rm = TRUE)
  # 
  # tmp=0
  # for (i in 1:length(Tm)){
  #   R=Tm[[i]]-offset[[i]]
  #   #Es=rbind(Es,matrix(as.numeric(R)))
  #   tmp=tmp+sum(as.numeric(R)^2,na.rm = TRUE)
  #   
  #   
  # }  
  
  N <- nrow(Tm)
  P <- ncol(Tm)

  if(is.null(nu1)){ nu1<-1 } 
  if(is.null(s10)){ s10<-1 } 
  
  1/rgamma(1, (N*P+s10)/2 , (tmp.1+nu1*s10)/2 )
}
