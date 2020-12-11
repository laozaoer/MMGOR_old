\name{MM_est}
\alias{MM_est}
\title{
Analysing clustered current status data}
\description{
MM_est is used to analyze clustered current status data. The function provides parameter estimates and corresponding variance estimates. 
}
\usage{
MM_est(initial_value,quad_mat,Delta,X,Z,n,ni,r,Design_mat,betadim,gammadim,itermax=500,tol=1e-7)

}
\arguments{
  \item{initial_value}{The initial value of the estimate}
  \item{quad_mat}{
            A matrix, Gauss-Hermite quadrature}
  \item{Delta}{
        A list, each component of the list is a vector containing  indicators of the censoring for each observation.}
  \item{X}{
  A list, each component is a matrix containing the within-subject level covariates.
  }
  \item{Z}{
  A matrix, contains the subject level covariates.
}
  \item{n}{Numbers of total subjects.}
  \item{ni}{
   Numbers of observations within each subject.}
  \item{r}{
  Specific r value for the generalized odds ratio model.}
  \item{Design_mat}{
  A list, each component is a matrix containing the spline function values for each observation.
    }
  \item{betadim}{
  Dimensions of within-subject level covariates.
    }
  \item{gammadim}{
  Dimensions of subject level covariates.
 .}
  \item{itermax}{
  Maximum iteration times, the default value is 500. 
  }
   \item{tol}{
   This denotes the summation of the absolute values of the relative tolerance of all parameters in the model. It is used to define the convergence of the parameter estimates. The default value is 1e-7.
   }
}

\value{
Function MM_est returns matrix. The first column contains the estimate of the regression parameters. The second column contains the variance estimates of these estimators.


}

