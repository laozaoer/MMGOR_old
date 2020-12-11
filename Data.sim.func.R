H=function(t){
  return(log(1+t)+t^(3/2))
}
Data_object=function(t,X,Z,b,r,beta,gamma,theta,U,H){
  Ht=H(t)
  if(r>0){
    S=(1+r*Ht*exp(beta*X+gamma*Z+theta*b))^(-1/r)
  }
  else{
    S=exp(-Ht*exp(beta*X+gamma*Z+theta*b))
  }
  
  return(S-U)
}
Generate_singleT=function(X,Z,b,r,beta,gamma,theta,H){
  U=runif(1,0,1)
  tresult=nleqslv::nleqslv(x=0.1,fn=Data_object,X=X,Z=Z,b=b,r=r,beta=beta,gamma=gamma,theta=theta,U=U,H=H)
  return(tresult$x)
}
Generate_T=function(X,Z,b,r,beta,gamma,theta,n,ni,H){
  result=list()
  length(result)=n
  for (i in 1:n) {
    result[[i]]=rep(0,ni[i])
    for (j in 1:ni[i]) {
      result[[i]][j]=Generate_singleT(X[[i]][j],Z[i],b[i],r,beta,gamma,theta,H)
    }
  }
  return(result)
}
data_for_est=function(r,beta,gamma,theta,X,Z,n,ni,knotsnum=2,order=2,H){
  rawC=Generate_T(X,Z,b,r,beta,gamma,theta,n,mi,H)
  lowC=0
  upC=quantile(as.numeric(as.character(unlist(rawC))),probs = 0.85)
  
  for(i in 1:n){
    C[[i]]=runif(mi[i],lowC,upC)
  }
  Delta=list()
  length(Delta)=n
  for (i in 1:n) {
    Delta[[i]]=rep(0,mi[i])
    for (j in 1:mi[i]) {
      if(rawC[[i]][j]<=C[[i]][j]){
        Delta[[i]][j]=1
      }
    }
  }
  
  blC <- list()
  length(blC) <- n
  knots <- seq(0,1  , length.out = (knotsnum + 2))
  knots=knots[3:length(knots)-1]
  for (i in 1:n) {
    blC[[i]]=t(ibs((C[[i]]-lowC)/(upC-lowC),knots = knots,degree=order,Boundary.knots = c(0,1),intercept = TRUE))
  }
  return(list(Delta=Delta,spline_value=blC))
}

library(gaussquad)
library(numDeriv)
library(ucminf)
library(extraDistr)
library(fda)
library(splines2)
library(nleqslv)