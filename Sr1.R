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
MM_est=function(initial_value,quad_mat,Delta,X,Z,n,ni,r,Design_mat,betadim,gammadim,itermax=500,tol=1e-7){
  parest=MainFunc(initial_value,quad_mat,Delta,X,Z,n,ni,r,Design_mat,betadim,gammadim,itermax,tol)
  hessian=numDeriv::hessian(func=testquadrature1current,x=parest[,1],rules=quad_mat,Delta=Delta,
                            X=X,Z=(Z),n=n,ni=ni,r=r,blC=Design_mat,betadim=betadim,gammadim=gammadim)
  var=diag(solve(hessian,tol=1e-40))
  result=cbind(parest[(1:(betadim+gammadim+1)),1],var[1:(betadim+gammadim+1)])
  result[(betadim+gammadim+1),1]=exp(result[(betadim+gammadim+1),1])
  result[(betadim+gammadim+1),2]=(result[(betadim+gammadim+1),1])^2*result[(betadim+gammadim+1),2]
  result=as.data.frame(result)
  colnames(result)=c("par.est","var.est")
  return(result)
}

library(gaussquad)
library(numDeriv)
library(ucminf)
library(extraDistr)
library(fda)
library(splines2)
library(nleqslv)
library(MMGOR)
r=0
set.seed(30000+2*19)
xi=0
beta=as.vector(-1)
gamma=as.vector(-1)
theta=1

n=300
# Z=rnorm(n,0,1)
Z=runif(n,-1,1)
mi=rep(0,n)
b=rnorm(n,0,1)
for(i in 1:n){
  mi[i]=rtpois(1,exp(1.7),a=1,b=8)
}
C=list()
length(C)=n
for(i in 1:n){
  # C[[i]]=rexp(n=mi[i],1)
  C[[i]]=runif(mi[i],0,1)
}
X=list()
length(X)=n
for (i in 1:n) {
  # X[[i]]=as.matrix(rbern(mi[i],0.578))
  X[[i]]=as.matrix(runif(mi[i],-1,1))
}

rawC=Generate_T(X,Z,b,r,beta,gamma,theta,n,mi)
lowC=0
upC=quantile(as.numeric(as.character(unlist(rawC))),probs = 0.85)

for(i in 1:n){
  # C[[i]]=rexp(n=mi[i],1)
  C[[i]]=runif(mi[i],lowC,upC)
}
Delta=list()
length(Delta)=n
for (i in 1:n) {
  Delta[[i]]=rep(0,mi[i])
  for (j in 1:mi[i]) {
    if(rawC[[i]][j]<=C[[i]][j]){
      Delta[[i]][j]=1
      # C[[i]][j]=rawC[[i]][j]
    }
  }
}
knotsnum=2
order=2
blC <- list()
length(blC) <- n
knots <- seq(0,1  , length.out = (knotsnum + 2))
knots=knots[3:length(knots)-1]
for (i in 1:n) {
  blC[[i]]=t(ibs((C[[i]]-lowC)/(upC-lowC),knots = knots,degree=order,Boundary.knots = c(0,1),intercept = TRUE))
}
Z=as.matrix(Z)

myrules=hermite.h.quadrature.rules(30,normalized=FALSE)
myrules=as.matrix(myrules[[30]])

# start_time <- Sys.time()
# result=MainFunc(rep(0,8),myrules,Delta,X,Z,n,mi,r,blC,1,1,500,1e-7)
# end_time <- Sys.time()

MM_est=function(initial_value,quad_mat,Delta,X,Z,n,ni,r,Design_mat,betadim,gammadim,itermax=500,tol=1e-7){
  parest=MainFunc(initial_value,quad_mat,Delta,X,Z,n,ni,r,Design_mat,betadim,gammadim,itermax,tol)
  hessian=numDeriv::hessian(func=testquadrature1current,x=parest[,1],rules=quad_mat,Delta=Delta,
                            X=X,Z=(Z),n=n,ni=ni,r=r,blC=Design_mat,betadim=betadim,gammadim=gammadim)
  var=diag(solve(hessian,tol=1e-40))
  result=cbind(parest[(1:(betadim+gammadim+1)),1],var[1:(betadim+gammadim+1)])
  result[(betadim+gammadim+1),1]=exp(result[(betadim+gammadim+1),1])
  result[(betadim+gammadim+1),2]=(result[(betadim+gammadim+1),1])^2*result[(betadim+gammadim+1),2]
  result=as.data.frame(result)
  colnames(result)=c("par.est","var.est")
  return(result)
}

MM_est(rep(0,8),myrules,Delta,X,Z,n,mi,r,blC,1,1)
