% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runStepwiseReg.R
\name{runStepwiseReg}
\alias{runStepwiseReg}
\title{Run step-wise regression to order the features}
\usage{
runStepwiseReg(ggc, filt.data)
}
\arguments{
\item{ggc}{gene-gene correlation matrix}

\item{filt.data}{filtered and normalised log-transformed genes x cells single-cell RNA-seq data matrix}
}
\value{
optimal feature set
}
\description{
Run step-wise regression to order the features
}
