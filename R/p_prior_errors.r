#' Computes prior probability of errors
#'
#' This is an auxiliary function in single package. It takes a data frame with counts by position, nucleotide and Qscore and it summarises it into proportion of nucleotide counts by position.
#' @param counts_pnq Data frame with columns position nucleoide quality counts, as returned by parse_countspnq
#' @param output_file File name for output, if save=TRUE.
#' @param save Logical. Should data be saved in a output_file?
#' @return Data frame with columns position nucleotide prior.error.
#' @importFrom utils write.table
#' @importFrom rlang .data
#' @import dplyr
#' @export p_prior_errors
#' @examples
#' refseq_fasta <- system.file("extdata", "ref_seq.fasta", package = "single")
#' train_reads_example <- system.file("extdata", "train_seqs_500.sorted.bam",
#'                                    package = "single")
#' counts_pnq <- pileup_by_QUAL(bam_file=train_reads_example,
#'     pos_start=1,pos_end=10)
#' p_prior_errors <- p_prior_errors(counts_pnq=counts_pnq)
#' head(p_prior_errors)
p_prior_errors         <- function(counts_pnq,output_file=NULL, save=FALSE){
    prior_error <- counts_pnq %>%
        dplyr::group_by(.data$pos,.data$strand, .data$nucleotide)%>%
        dplyr::summarise(count=sum(.data$count))%>%
        dplyr::ungroup() %>%
        dplyr::group_by(.data$strand,.data$pos)%>%
        dplyr::mutate(p_prior_error = .data$count/sum(.data$count)) %>%
        dplyr::ungroup()%>%
        dplyr::select(.data$strand,.data$pos, .data$nucleotide, .data$p_prior_error)

    if(save){utils::write.table(prior_error, file = output_file, row.names = FALSE)}
    return(prior_error)
}
