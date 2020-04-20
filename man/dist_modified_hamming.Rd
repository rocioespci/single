\name{dist_modified_hamming}
\alias{dist_modified_hamming}
\title{dist_modified_hamming}
\usage{
dist_modified_hamming(seq_1,seq_2,ref, w_1, w_2)
}
\description{
Distance between strings. Each input must have the same length.
}
\arguments{
  \item{seq_1}{character vector}
  \item{seq_2}{character vector}
  \item{ref}{character vector}
  \item{w_1}{numeric vector}
  \item{w_1}{numeric vector}
}
\value{
Numeric, distance between seq_1 and seq_2, weighted by w_1 and w_2 and in comparisson to reference sequence ref.

For each position in the vector, it sums:
\itemize{
\item if seq_1 and seq_2 are equal to reference, 0
\item if seq_1 and seq_2 are different to reference and different between them, -w_1  w_2
\item if seq_1 and seq_2 are different to reference and equal between them, +w_1 w_2
}


}
