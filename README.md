## Overview
The **MMGOR** package implements a novel minorize-maximize algorithm to estimate parameters of generalized odds rate (GOR) model with clustered current status data. **MMGOR** allows any nonnegative r values of GOR model, which covers a wide range of commonly used survival models. The package takes advantage of C++ computational efficiency to reduce computation time.

## Installation

Downlowd the compressed file *MMGOR_1.0.tar.gz* to the working directory and  use following command to install it
```
install.packages("MMGOR_1.0.tar.gz")
```

## Example
There is an example to show how to use this package. First, use following command to simulate data. Data simulating has several dependencies [**gaussquad**](https://cran.r-project.org/web/packages/gaussquad/index.html), [**numDeriv**](https://cran.r-project.org/web/packages/numDeriv/index.html), [**extraDistr**](https://cran.r-project.org/web/packages/extraDistr/index.html), [**fda**](https://cran.r-project.org/web/packages/fda/index.html), [**splines2**](https://cran.r-project.org/web/packages/splines2/index.html) and [**nleqslv**](https://cran.r-project.org/web/packages/nleqslv/index.html).
```
# Load package
library(MMGOR)
```
### Scenario 1: <img src="http://chart.googleapis.com/chart?cht=tx&chl= r=0" style="border:none;">, <img src="http://chart.googleapis.com/chart?cht=tx&chl= \theta=1" style="border:none;">, <img src="http://chart.googleapis.com/chart?cht=tx&chl= \beta=-1" style="border:none;">, <img src="http://chart.googleapis.com/chart?cht=tx&chl= \gamma=-1" style="border:none;">, <img src="http://chart.googleapis.com/chart?cht=tx&chl= n=300" style="border:none;">

#### Simulate data:
```
library(MMGOR)
set.seed(718)
H=function(t) log(1+t)+t^(3/2) 
data=data_for_est(r=0,beta=c(-1),gamma=c(-1),theta=1,n=300,H=H,knotsnum=2,order=2,quadnum=30)
```
This function generates 300 subjects, where each subject has up to 8 observations. Both the covariates of subject level and within-subject level are generated from uniform(-1,1) distribution. Users can use their own H function by changing the specific H function form.

#### Estimate parameters:
Regression parameters can be estimated using the command
```
result=MM_est(rep(0,8),data$GHrules,data$Delta,data$X,data$Z,data$n,data$ni,data$r,data$spline_value,data$betadim,data$gammadim)
```
The output includes the parameters estimates and the corresponding variance estimates.
