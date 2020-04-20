\name{modified_hamming_distance_matrix}
\alias{modified_hamming_distance_matrix}
\title{modified_hamming_distance_matrix}
\usage{
modified_hamming_distance_matrix(filename_to_evaluate, reference_sequence,weighting=c("raw","naive","priors"),dist_function="dist_modified_hamming", cut_off=0.9,filename_out=NULL,filename_ids=NULL)
}
\description{
HIgh level function to compute distance matrix between several sequences, stored in individual files in each filename_to_evaluate path.}
\arguments{
  \item{filename_to_evaluate}{character vector, paths where individual files are stored.}
  \item{reference_sequence}{character vector, with reference sequence.}
  \item{weighting}{Which weights should be used for distance calculations. One of "raw","naive","priors".}
  \item{dist_function}{Character, distance function to be used. One of "dist_modified_hamming","dist_different_minus_equal", "dist_modified_hamming_soft","dist_equal".}
  \item{cut_off}{numeric between 0 and 1 (length 1).}
  \item{filename_out}{Character. Name of file in which distacne matrix should be stored. If NULL, it will not be saved.}
  \item{filename_ids}{If length(filename_to_evaluate)>1, a vector of length(filename_to_evaluate) giving individual labels to each filename_to_evaluate.}
}
\value{
Numeric distance matrix between sequences and in comparisson to reference sequence.

}
