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

# Set the true value of parameters
r=0 
beta=as.vector(-1)
gamma=as.vector(-1)
theta=1
n=300

# Generate subject level covariates. Users can change the distribution of covariates.
Z=runif(n,-1,1)

# Generate numbers of observations within each subjects
mi=rep(0,n) 
for(i in 1:n){
  mi[i]=rtpois(1,exp(1.7),a=1,b=8)
}

# Generate cluster effect
b=rnorm(n,0,1)
C=list()
length(C)=n
for(i in 1:n){
  C[[i]]=rep(0,mi[i])
}

# Generate within-subject level covariates
X=list()
length(X)=n
for (i in 1:n) {
  X[[i]]=as.matrix(runif(mi[i],-1,1))
}
Z=as.matrix(Z)

# Generate Gauss-Hermite Quadrature rule
myrules=hermite.h.quadrature.rules(30,normalized=FALSE)
myrules=as.matrix(myrules[[30]])

# Generate the censoring indicators and the spline function values. Users can use their own H function.
data=data_for_est(r,beta,gamma,theta,X,Z,n,mi,knotsnum=2,order=2,H)

# Use the proposed MM algorithm to obtain the estimates.

result=MM_est(rep(0,8),myrules,data[[1]],X,Z,n,mi,r,data[[2]],1,1)

```

### Scenario 2: <img src="http://chart.googleapis.com/chart?cht=tx&chl= r=2" style="border:none;">, <img src="http://chart.googleapis.com/chart?cht=tx&chl= \theta=1" style="border:none;">, <img src="http://chart.googleapis.com/chart?cht=tx&chl= \beta=-1" style="border:none;">, <img src="http://chart.googleapis.com/chart?cht=tx&chl= \gamma=-1" style="border:none;">, <img src="http://chart.googleapis.com/chart?cht=tx&chl= n=300" style="border:none;">
```
# Simulate data 

# Set the true value of parameters
r=2
beta=as.vector(-1)
gamma=as.vector(-1)
theta=1
n=300

# Generate subject level covariates
Z=runif(n,-1,1)

# Generate numbers of observations within each subjects
mi=rep(0,n) 
for(i in 1:n){
  mi[i]=rtpois(1,exp(1.7),a=1,b=8)
}

# Generate cluster effect
b=rnorm(n,0,1)
C=list()
length(C)=n
for(i in 1:n){
  C[[i]]=rep(0,mi[i])
}

# Generate within-subject level covariates
X=list()
length(X)=n
for (i in 1:n) {
  X[[i]]=as.matrix(runif(mi[i],-1,1))
}
Z=as.matrix(Z)

# Generate Gauss-Hermite Quadrature rule
myrules=hermite.h.quadrature.rules(30,normalized=FALSE)
myrules=as.matrix(myrules[[30]])

# Generate the censoring indicators and the spline function values. Users can use different H function.
data=data_for_est(r,beta,gamma,theta,X,Z,n,mi,knotsnum=2,order=2,H)

# Use the proposed MM algorithm to obtain the estimates.

result=MM_est(rep(0,8),myrules,data[[1]],X,Z,n,mi,r,data[[2]],1,1)

```

