% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/testCCF.R
\name{calCCF}
\alias{calCCF}
\title{Calulate the CCF second test value}
\usage{
calCCF(ccfResult)
}
\arguments{
\item{ccfResult}{The result of CCF function}
}
\value{
A number of CCF sencond test value
}
\description{
Takes in the result of CCF
}
\details{
CCF second value = (max(cffResult(ro>0))-max(cffResult(ro<0)))/(max(cffResult(ro>0))+max(cffResult(ro<0)))
}

