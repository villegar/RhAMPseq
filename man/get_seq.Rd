% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_seq}
\alias{get_seq}
\title{Combine list of bases}
\usage{
get_seq(raw, upper = TRUE)
}
\arguments{
\item{raw}{list of bases}

\item{upper}{flag on whether the sequence should be returned in all caps}
}
\value{
string of bases
}
\description{
Combine list of bases
}
\examples{
set.seed(123)
# Generate random list of DNA bases
seq_list <- list(c("A", "C", "G", "T")[sample(1:4, 100, TRUE)])
get_seq(seq_list) 
get_seq(seq_list, FALSE) # Return lower case sequence
}
