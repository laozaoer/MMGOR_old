## MMGOR
A novel minorize-maximize algorithm for the generalized odds ratio Model for clustered current status data

## Installation

To use the proposed MM algorithm, users should first downlowd the compressed file *MMGOR_1.0.tar.gz* and install it.
```
install.packages("MMGOR_1.0.tar.gz")
```

## Examples
```
# Load package
library(MMGOR)

#Import required packages and functions for generating data
source("Data.sim.func.R")
```
### Scenario 1: <img src="http://chart.googleapis.com/chart?cht=tx&chl= r=0" style="border:none;">, <img src="http://chart.googleapis.com/chart?cht=tx&chl= \theta=1" style="border:none;">, <img src="http://chart.googleapis.com/chart?cht=tx&chl= \beta=-1" style="border:none;">, <img src="http://chart.googleapis.com/chart?cht=tx&chl= \gamma=-1" style="border:none;">, <img src="http://chart.googleapis.com/chart?cht=tx&chl= n=300" style="border:none;">
```
# Simulate data 
r=0 
beta=as.vector(-1)
gamma=as.vector(-1)
theta=1
n=300
# generating subject level covariates
Z=runif(n,-1,1)
mi=rep(0,n)// numbers of observations within each subjects
for(i in 1:n){
  mi[i]=rtpois(1,exp(1.7),a=1,b=8)
}
b=rnorm(n,0,1)// simulated clustered effect
C=list()
length(C)=n
for(i in 1:n){
  C[[i]]=rep(0,mi[i])
}
X=list()// within-subject level covariates
length(X)=n
for (i in 1:n) {
  X[[i]]=as.matrix(runif(mi[i],-1,1))
}
Z=as.matrix(Z)

myrules=hermite.h.quadrature.rules(30,normalized=FALSE)
myrules=as.matrix(myrules[[30]])
data=data_for_est(r,beta,gamma,theta,X,Z,n,mi,knotsnum=2,order=2,H)
```

### Simulate data with $r=1$, $\theta=1$, $\beta=-1$, $\gamma=-1$, $n=300$  and estimate the parameters.
```
source("Sr1.R")
result=MM_est(rep(0,8),myrules,Delta,X,Z,n,mi,r,blC,1,1)
```

### Simulate data with $r=2$, $\theta=1$, $\beta=-1$, $\gamma=-1$, $n=300$  and estimate the parameters.
```
source("Sr2.R")
result=MM_est(rep(0,8),myrules,Delta,X,Z,n,mi,r,blC,1,1)
```
