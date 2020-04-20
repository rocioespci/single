\name{accomodate_dist}
\alias{accomodate_dist}
\title{accomodate_dist}
\usage{
accomodate_dist(x, names=NULL)
}
\description{
Builds a symmetric, named matrix.
}
\arguments{
  \item{x}{Numeric squared matrix.}
  \item{names}{Character vector, length equal to number of columns of x.}
}
\value{
Matrix equal to x+t(x) with rows and cols names equal to names.

}
