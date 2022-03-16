#' Computes prior probability of mutations
#'
#' This is an auxiliary function in single package. It computes the prior probability of mutation in a gene library.
#' @param rates.matrix Mutation rate matrix: 4x5 matrix, each row/col representing a nucleotide (col adds deletion), and the values is the mutational rate from row to col.
#' @param mean.n.mut Mean number of mutations expected (one number).
#' @param ref_seq DNAStringSet containing the true reference sequence.
#' @param save Logical. Should data be saved in a output_file?
#' @param output_file File name for output, if save=TRUE.
#' @return Data frame with columns wt.base (wild type nucleotide), nucleotide (mutated nucleotide), p_mutation (probaility of mutation)
#' @importFrom reshape2 melt
#' @importFrom utils write.table
#' @importFrom Biostrings oligonucleotideFrequency
#' @export p_prior_mutations
#' @examples
#' pos_start <- 1
#' pos_end <- 10
#' refseq_fasta <- system.file("extdata", "ref_seq.fasta", package = "single")
#' train_reads_example <- system.file("extdata", "train_seqs_500.sorted.bam",
#'                                    package = "single")
#' counts_pnq <- pileup_by_QUAL(bam_file=train_reads_example,
#'     pos_start=pos_start,pos_end=pos_end)
#' p_prior_mutations <- p_prior_mutations(rates.matrix = mutation_rate,
#'     mean.n.mut = 5,ref_seq = ref_seq,save = FALSE)
#' head(p_prior_mutations)

p_prior_mutations      <- function(rates.matrix, mean.n.mut, ref_seq,
                            save=FALSE, output_file="tablePriorMutations.txt"){
    composition_wt <- c(Biostrings::oligonucleotideFrequency(ref_seq,1),
                        length(ref_seq))
    names(composition_wt)[5] <- "-"

    mutations_rate              <- apply(rates.matrix, 2, sum, na.rm=TRUE)
    expected_mutations_perbase  <- composition_wt*mutations_rate
    normalization_factor        <- sum(expected_mutations_perbase)/mean.n.mut

    expected_mutation_rate <- rates.matrix / normalization_factor
    expected_mutation_rate <- reshape2::melt(expected_mutation_rate,
                                        varnames = c("wt.base","nucleotide"),
                                        value.name = "p_mutation")
    expected_mutation_rate <- expected_mutation_rate[!is.na(expected_mutation_rate$p_mutation),]
    if(save){
        utils::write.table(expected_mutation_rate,
                file = output_file, row.names = FALSE)
    }
    return(expected_mutation_rate)
}
