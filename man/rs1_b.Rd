% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rs1_b.R
\name{rs1_b}
\alias{rs1_b}
\title{Gibbs update for behavior variance}
\usage{
rs1_b(Tm, offset=0,nu1=NULL,s10=NULL)
}
\arguments{
\item{Tm}{a list of V X P normal behavior matrix}

\item{nu1}{prior degrees of freedom}

\item{s10}{prior estimate of s1}
}
\value{
a new value of s1
}
\description{
Gibbs update for behavior variance
}
\author{
Selena Wang
}
