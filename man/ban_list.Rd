% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funs_data.R
\docType{data}
\name{ban_list}
\alias{ban_list}
\title{Vector of words banned to appear in the name of a pathway}
\format{
Vector of words banned to appear in the name of a pathway
\describe{
  \item{ban_list}{Vector of words banned to appear in the name of a pathway}
}
}
\source{
{
  \url{experience}
}
}
\usage{
ban_list
}
\description{
GO, KEGG, MSIGDB include pathways which are not informative or very context specific that would be
strange to use in most of the analysis, this list is removed from the database
}
\keyword{datasets}
