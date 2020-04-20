\name{dist_modified_hamming_soft}
\alias{dist_modified_hamming_soft}
\title{dist_modified_hamming_soft}
\usage{
dist_modified_hamming_soft(seq_1,seq_2,ref, w_1, w_2, cut_off=0.9)
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
  \item{cut_off}{numeric (length 1)}
}
\value{
Numeric, distance between seq_1 and seq_2, in comparisson to reference sequence ref.

Calculated as:
\itemize{
\item If both seq_1 and seq_2 different to reference and between them, and w_1 or w_2 over cut_off, then sum w_1 w_2.
\item Else, sums 0.
}
}
\seealso{
\code{\link{dist_modified_hamming}}
}

