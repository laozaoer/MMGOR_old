Data_object=function(t,X,Z,b,r,beta,gamma,theta,U){
  Ht=log(1+t)+t^(3/2)
  if(r>0){
    S=(1+r*Ht*exp(beta*X+gamma*Z+theta*b))^(-1/r)
  }
  else{
    S=exp(-Ht*exp(beta*X+gamma*Z+theta*b))
  }
  
  return(S-U)
}
Generate_singleT=function(X,Z,b,r,beta,gamma,theta){
  U=runif(1,0,1)
  tresult=nleqslv::nleqslv(x=0.1,fn=Data_object,X=X,Z=Z,b=b,r=r,beta=beta,gamma=gamma,theta=theta,U=U)
  return(tresult$x)
}
Generate_T=function(X,Z,b,r,beta,gamma,theta,n,ni){
  result=list()
  length(result)=n
  for (i in 1:n) {
    result[[i]]=rep(0,ni[i])
    for (j in 1:ni[i]) {
      result[[i]][j]=Generate_singleT(X[[i]][j],Z[i],b[i],r,beta,gamma,theta)
    }
  }
  return(result)
}
library(gaussquad)
library(numDeriv)
library(ucminf)
library(extraDistr)
library(fda)
library(splines2)
library(nleqslv)