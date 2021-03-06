% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spacemut_functions.R
\name{select.ncomponents}
\alias{select.ncomponents}
\title{Estimate number of components with reflected property and processes for different input extracted components}
\usage{
select.ncomponents(rate, n.min = 2, n.max = 30, cutoff = 0.8,
  n.cores = 1)
}
\arguments{
\item{rate}{matrix of mutation rates of mutation types in genomic windows}

\item{n.min}{min number of components to extract}

\item{n.max}{max number of components to extract}

\item{cutoff}{cutoff on reflection property}

\item{n.cores}{nunmber of cores to use}
}
\value{
number of components (with reflection property) and processes per each input number of extracted components
}
\description{
Estimate number of components with reflected property and processes for different input extracted components
}
