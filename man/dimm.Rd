% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimation-funcs-package.R
\name{dimm}
\alias{dimm}
\title{fit a regression model with high dimensional correlated outcome}
\usage{
dimm(formula, data, id = id, response_indicator = NULL,
  subject_indicator = NULL, family, corstr, method = NULL, lam = NULL,
  folds = NULL, cluster = NULL, ...)
}
\arguments{
\item{formula}{an object of class "formula" of the type \code{Y ~ X} indicating the response ~ covariate relationship to model}

\item{data}{a data frame containing the variables in the model}

\item{id}{subject id of participants}

\item{response_indicator}{a vector of integers from 1 up to M, indicating which block the response belongs to}

\item{subject_indicator}{a vector of integers from 1 up to N, indicating which subject block the subject belongs to. Default is one block}

\item{family}{"normal" for continuous response}

\item{corstr}{correlation structure to model the covariance of the response within each block}

\item{method}{method of combination to be used: "iterative" for iterative minimization of the objective function; "exact" for one-step estimator; "ridge" for regularized weight matrix}

\item{lam}{optional value of tuning parameter for regularization when method="ridge"; if supplied a vector, selection of best tuning parameter done by cross-validation}

\item{folds}{optional number of folds to be used in cross-validation selection of the tuning parameter when method="ridge"}

\item{cluster}{optional number of threads for parallelization of block analyses}

\item{...}{additional arguments to be passed}
}
\value{
an object of class "dimm". The function summary is used to obtain and print a summary of the results. The generic accessor functions coefficients, vcov and fitted.values extract features of the object.
}
\description{
Divides the data into blocks, analyzes blocks using composite likelihood, and combines results using GMM
}
\details{
A formula statement has an implicit intercept. To remove the intercept, use either \code{Y ~ X-1} or \code{X ~ 0+X}. dimm calls functions dimm.compute.mean or dimm.compute.mean.parallel for actual computation.
}
