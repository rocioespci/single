\name{clusters_consensus}
\alias{clusters_consensus}
\title{clusters_consensus}
\usage{
clusters_consensus(sequences,clusters_vector,reference_sequence,weights_vector=NULL,cutoff_prob=0.2)
}
\description{
Computes a consensus sequence for each cluster.
}
\arguments{
  \item{sequences}{Character matrix. Sequences over which compute consensus. Each row is one sequcene.}
  \item{clusters_vector}{Numeric vector. Indicates to which cluster belongs each sequence. As returned by cutree. It must have the same order as sequences' rows.}
  \item{reference_sequence}{Character vector. Reference sequence.}
  \item{weights_vector}{Numeric matrix. Probabily of being correct for each nucleotide in sequences, in the same shape and order. If NULL all be equal to one.}
  \item{cutoff_prob}{Numeric. A cut-off for which nucleotides in sequences will be ignored if weight_vector is below it.}
}
\value{
Returns a two columns data frame: Cluster number and mutations respect to reference_sequence in the consensus sequence.
}
