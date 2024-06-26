% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multinom_model.R
\name{fit_multinom_model}
\alias{fit_multinom_model}
\title{Fit a multinomial modelto predict annotation using NMF factors as covariates.}
\usage{
fit_multinom_model(factors, source_annotations)
}
\arguments{
\item{factors}{a matrix of factors from NMF}

\item{source_annotations}{a vector of annotations}
}
\value{
A multinomial model object
}
\description{
Fit a multinomial modelto predict annotation using NMF factors as covariates.
}
