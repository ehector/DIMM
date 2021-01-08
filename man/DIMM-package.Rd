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
