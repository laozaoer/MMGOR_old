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
data_for_est=function(r,beta,gamma,theta,n,H,knotsnum=2,order=2,quadnum=30){
  
  Z=runif(n,-1,1)
  mi=rep(0,n)
  b=rnorm(n,0,1)
  for(i in 1:n){
    mi[i]=rtpois(1,exp(1.7),a=1,b=8)
  }
  C=list()
  length(C)=n
  for(i in 1:n){
    C[[i]]=runif(mi[i],0,1)
  }
  X=list()
  length(X)=n
  for (i in 1:n) {
    X[[i]]=as.matrix(runif(mi[i],-1,1))
  }
  Z=as.matrix(Z)
  
  myrules=hermite.h.quadrature.rules(quadnum,normalized=FALSE)
  myrules=as.matrix(myrules[[quadnum]])
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
  return(list(X=X,Z=Z,n=n,ni=mi,r=r,Delta=Delta,spline_value=blC,GHrules=myrules,betadim=length(beta),gammadim=length(gamma)))
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
