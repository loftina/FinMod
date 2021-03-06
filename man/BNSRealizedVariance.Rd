% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BNS.R
\name{BNSRealizedVariance}
\alias{BNSRealizedVariance}
\title{Barndorff‐Nielsen and Shephard Realized Variance}
\usage{
BNSRealizedVariance(realizedVar_0, time, lambda, delta, gamma, phi)
}
\arguments{
\item{realizedVar_0}{Initial realized variance of a stock repeated T times}

\item{time}{Number of realized variance observations}

\item{lambda}{The mean reverting variable}

\item{delta}{Influences the tail thickness of the disribution}

\item{gamma}{Influences the scale of the distribution}

\item{phi}{Volatility of the volatility.}
}
\value{
Realized variance according to the Barndorff‐Nielsen and Shephard model when the correct kappa value is used
}
\description{
The Barndorff‐Nielsen and Shephard model's approximation for realized variance.
}
