% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funs_similarities.R
\name{get_AUtri}
\alias{get_AUtri}
\title{Triangle area between three points}
\usage{
get_AUtri(b_value, c_value)
}
\arguments{
\item{b_value}{Value}

\item{c_value}{Value}
}
\value{
Triangle area
}
\description{
One is fixed (0,0), the second is (b_value,0), the third is (0,c_value)
Triangle area increases more when two points increase
Triangle area penalizes when only one point changes
}
