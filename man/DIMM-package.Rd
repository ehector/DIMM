\name{DIMM-package}
\alias{DIMM-package}
\alias{DIMM}
\docType{package}
\title{
\packageTitle{DIMM}
}
\description{
\packageDescription{DIMM}
}
\details{
  
  THIS SOURCE CODE IS SUPPLIED "AS IS" WITHOUT WARRANTY OF ANY KIND, AND ITS AUTHOR DISCLAIMS ANY AND ALL WARRANTIES, INCLUDING BUT NOT LIMITED TO ANY IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, AND ANY WARRANTIES OR NON INFRINGEMENT. THE USER ASSUMES ALL LIABILITY AND RESPONSIBILITY FOR USE OF THIS SOURCE CODE, AND THE AUTHOR WILL NOT BE LIABLE FOR DAMAGES OF ANY KIND RESULTING FROM ITS USE. Without limiting the generality of the foregoing, the author does not warrant that the Source Code will be error-free, will operate without interruption, or will meet the needs of the user.

The DESCRIPTION file:
\packageDescription{DIMM}
\packageIndices{DIMM}
This easy to use package runs the Distributed and Integrated Method of Moments (DIMM) for regression with high dimensional outcomes. The user supplies a model and indicator vector. DIMM divides the data into blocks according to the supplied indicator, analyzes blocks using composite likelihood, and combines results using an optimal generalized method of moments (GMM) equation.
}
\author{
\packageAuthor{DIMM}

Maintainer: \packageMaintainer{DIMM}
}
\references{
The posdef.matrix function was written by Ravi Varadhan: https://stat.ethz.ch/pipermail/r-help/2008-February/153708.
  
Emily C. Hector and Peter X.-K. Song (2020). A distributed and integrated method of moments for high-dimensional correlated data analysis. Journal of the American Statistical Association, pages 1–14. doi: 10.1080/01621459.2020.1736082.

Lars P. Hansen. (1982), Large sample properties of generalized method of moments estimators, Econometrica 50(4), 1029-1054. \cr

Christiano Varin, Nancy Reid & David Firth. (2011), An overview of composite likelihood methods, Statistica Sinica 21(1), 5-42. \cr
}

\keyword{ package }

\examples{
########### Generate Data ########### 
N <- 1000                 # sample size
M <- 200                  # response dimension
family <- "normal"        # response is Gaussian distributed
method <- "exact"         # find one-step estimator
theta <- c(0.3, 0.6)      # the true parameter value

# The response indicator partitions the response into blocks
# The subject indicator partitions subjects into blocks
response_indicator <- c(rep(1,45), rep(2,42), rep(3,50), rep(4,34), rep(5,29))
subject_indicator <- rep(1,N)

# Generate a dataset with one covariance, response of 
# dimension M with block AR-1 covariance structure
data <- dataset.normal.X1(family, N, M, theta, response_indicator, seed=500)

########### Analysis with DIMM ###########

# Model block covariance with independence structure
results.ind <- dimm(response~X1, data, id="id", response_indicator, subject_indicator, 
                    family, corstr="independence", method, cluster=5)

# Model block covariance with compound symmetry structure
results.CS <- dimm(response~X1, data, id="id", response_indicator, subject_indicator, 
                   family, corstr="CS", method, cluster=5)

# Model block covariance with AR-1 structure                             
results.AR1 <- dimm(response~X1, data, id="id", response_indicator, subject_indicator, 
                    family, corstr="AR1", method, cluster=5)

# Model block covariance with AR-1 structure, 
# using regularized weight matrix                            
results.AR1 <- dimm(response~X1, data, id="id", response_indicator, subject_indicator, 
                    family, corstr="AR1", method="ridge", 
                    lam=seq(0,1,0.1), folds=5, cluster=5)
                                        
## Goodness-of-fit test
objective.dimm(results.AR1)

}
