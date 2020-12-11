## MMGOR
Minorize-Maximize Algorithm for the Generalized Odds Ratio Model for Clustered Current Status Data

## Installation

To use the proposed MM algorithm, users should first downlowd the compressed file *MMGOR_1.0.tar.gz* and install it.
```
install.packages("MMGOR_1.0.tar.gz")
```

## Examples
```
# Load package
library(MMGOR)

#Import required packages and functions
source("Data.sim.func.R")
```
### Simulate data with $r=0$, $\theta=1$, $\beta=-1$, $\gamma=-1$, $n=300$  and estimate the parameters.
```
source("Sr0.R")
result=MM_est(rep(0,8),myrules,Delta,X,Z,n,mi,r,blC,1,1)
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
