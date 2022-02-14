#' Computes prior probability of mutations
#'
#' This is an auxiliary function in single package. It computes the prior probability of mutation in a gene library.
#' @param rates.matrix Mutation rate matrix: 4x5 matrix, each row/col representing a nucleotide (col adds deletion), and the values is the mutational rate from row to col.
#' @param mean.n.mut Mean number of mutations expected (one number).
#' @param reference_sequence Reference sequence: vector of characters, as returned by load_reference_sequence
#' @param save Logical. Should data be saved in a output_file?
#' @param output_file File name for output, if save=TRUE.
#' @return Data frame with columns wt.base (wild type nucleotide), nucleotide (mutated nucleotide), p_mutation (probaility of mutation)
#' @importFrom reshape2 melt
#' @importFrom utils write.table
#' @export calculate_prior_mutations
#' @examples
#' ref_seq_file = system.file("extdata", "ref_seq.fasta", package = "single")
#' ref_seq = load_reference_sequence(ref_seq_file)
#' calculate_prior_mutations(mutation_rate, 3, ref_seq)
calculate_prior_mutations      <- function(rates.matrix, mean.n.mut, reference_sequence, save=FALSE, output_file="tablePriorMutations.txt"){
    composition_wt <- c(table(reference_sequence), length(reference_sequence))
    names(composition_wt)[5] <- "-"

    mutations_rate              <- apply(rates.matrix, 2, sum, na.rm=TRUE)
    expected_mutations_perbase  <- composition_wt*mutations_rate
    normalization_factor        <- sum(expected_mutations_perbase)/mean.n.mut

    expected_mutation_rate <- rates.matrix / normalization_factor
    expected_mutation_rate <- reshape2::melt(expected_mutation_rate, varnames = c("wt.base","nucleotide"), value.name = "p_mutation")
    expected_mutation_rate <- expected_mutation_rate[!is.na(expected_mutation_rate$p_mutation),]
    if(save){utils::write.table(expected_mutation_rate, file = output_file, row.names = FALSE)}
    return(expected_mutation_rate)
}
