## MMGOR
A novel minorize-maximize algorithm for the generalized odds ratio Model for clustered current status data

## Installation

Downlowd the compressed file *MMGOR_1.0.tar.gz* to the working directory and then install it.
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
H=function(t) log(1+t)+t^(3/2) // Users can use their own H function.
data=data_for_est(r=0,beta=c(-1),gamma=c(-1),theta=1,n=300,H=H,knotsnum=2,order=2,quadnum=30)
```
This function generates 300 subjects, where each subject has up to 8 observations. Both the covariates of subject level and within-subject level are generated from uniform(-1,1) distribution.

