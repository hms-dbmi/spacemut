% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spacemut_functions.R
\name{spectra.bootstrap}
\alias{spectra.bootstrap}
\title{Estimate robustness of components spectra using bootstrap}
\usage{
spectra.bootstrap(icM, rate, n.comp = ncol(icM), n.boot = 100,
  n.cores = 1)
}
\arguments{
\item{icM}{matrix of components spectra}

\item{rate}{matrix of mutation rates of mutation types in genomic windows}

\item{n.comp}{number of extract components in bootstrap. By default ,the same as in icM}

\item{n.boot}{number of bootstraps}

\item{n.cores}{number of cores to use}
}
\value{
matrix of Spearman correlations between original and bootstrapped spectra of components
}
\description{
Estimate robustness of components spectra using bootstrap
}
