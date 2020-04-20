\name{dist_different_minus_equal}
\alias{dist_different_minus_equal}
\title{dist_different_minus_equal}
\usage{
dist_different_minus_equal(seq_1,seq_2,ref)
}
\description{
Distance between strings. Each input must have the same length.
}
\arguments{
  \item{seq_1}{character vector}
  \item{seq_2}{character vector}
  \item{ref}{character vector}
}
\value{
Numeric, distance between seq_1 and seq_2, in comparisson to reference sequence ref.

Calculated as: number of positions in which seq_1 and seq_2 are different to reference and DIFERENT among them, minus number of positions in which seq_1 and seq_2 are different to reference and EQUAL among them.

}
