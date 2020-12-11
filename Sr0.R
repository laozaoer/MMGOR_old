r=0
beta=as.vector(-1)
gamma=as.vector(-1)
theta=1

n=300
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


