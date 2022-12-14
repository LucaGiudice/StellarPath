% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funs_counts.R
\name{get_rank_01}
\alias{get_rank_01}
\title{Apply 0-1 standardized rank transformation}
\usage{
get_rank_01(m, type = "total_rank")
}
\arguments{
\item{m}{Numeric matrix of count values (e.g. molecules x samples)}

\item{type}{Default sample, character string which defines the direction of ranking, it can be "sample_rank", "total_rank", "expression"}
}
\value{
0-1 standardized ranks of the counts for each sample (default) or for the total matrix (type: total_rank) or no transformation (type: expression)
}
\description{
This function transforms the sample's counts in 0-1 standardized ranks.
}
\details{
By default is type "sample_rank": For each sample, a molecule is ranked based on its count value and then the ranks are 0-1 standardized.
The most expressed molecule of a sample has the highest rank (e.g. 20000) divided by the number of molecules get 1
The least expressed molecule get a value equal to 0.

If type is "total_rank": a molecule is ranked based on its count value with respect all the count values in the matrix and then
the ranks are 0-1 standardized.

If type is "expression": no transformation is applied, lcpm count values are used
}
