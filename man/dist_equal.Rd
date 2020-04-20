\name{dist_equal}
\alias{dist_equal}
\title{dist_equal}
\usage{
dist_equal(seq_1,seq_2,ref, w_1, w_2, cut_off=0.9)
}
\description{
Distance between strings. Each input must have the same length (except cut_off).
}
\arguments{
  \item{seq_1}{character vector}
  \item{seq_2}{character vector}
  \item{ref}{character vector}
  \item{w_1}{numeric vector}
  \item{w_1}{numeric vector}
  \item{cut_off}{numeric between 0 and 1 (length 1).}
}
\value{
Numeric, distance between seq_1 and seq_2, weighted by w_1 and w_2 and in comparisson to reference sequence ref.

Number of positions in which seq_1, seq_2 and reference are all different and w_1 and w_2 are larger than cut_off.

}
