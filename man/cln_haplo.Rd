% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{cln_haplo}
\alias{cln_haplo}
\title{Clean haplotypes}
\usage{
cln_haplo(hap, read_length = 5, missing = NA)
}
\arguments{
\item{hap}{haplotype data with the format: \code{X/Y:n,m}.}

\item{read_length}{minimum read length to keep haplotype.}

\item{missing}{default value for haplotypes below the \code{read_length}.}
}
\value{
haplotype
}
\description{
Search for haplotypes below a threshold for read length and set them to
\code{missing}.
}
\examples{
cln_haplo("./.:0")
cln_haplo("4/2:34,24")
cln_haplo("2/2:71")
cln_haplo("4/2:2,1")
cln_haplo("2/2:4")
}
