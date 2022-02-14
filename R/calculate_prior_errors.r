#' Computes prior probability of errors
#'
#' This is an auxiliary function in single package. It takes a data frame with counts by position, nucleotide and Qscore and it summarises it into proportion of nucleotide counts by position.
#' @param counts_pnq Data frame with columns position nucleoide quality counts, as returned by parse_nucleotides_per_qscore.
#' @param output_file File name for output, if save=TRUE.
#' @param reference_sequence Reference sequence: vector of characters, as returned by load_reference_sequence
#' @param save Logical. Should data be saved in a output_file?
#' @return Data frame with columns position nucleotide prior.error.
#' @importFrom utils write.table
#' @importFrom rlang .data
#' @export calculate_prior_errors
#' @examples
#' ref_seq_file = system.file("extdata", "ref_seq.fasta", package = "single")
#' ref_seq = load_reference_sequence(ref_seq_file)
#' counts_pnq_file = system.file("extdata", "example_train_countsPNQ_parsed.txt", package = "single")
#' counts_pnq = read.table(file=counts_pnq_file, header=TRUE)
#' calculate_prior_errors(counts_pnq,reference_sequence=ref_seq)
calculate_prior_errors         <- function(counts_pnq,output_file=NULL, reference_sequence,save=FALSE){
    if(save){if(is.null(output_file)){stop('single::calculate_prior_errors Invalid output_file!')}}
    col_names <- c("position","nucleotide","quality","counts")
    if(any(colnames(counts_pnq)!=col_names)){
        stop("calculate_prior_errors: expected columns for counts_pnq are position nucleotide quality counts")
    }
    colnames(counts_pnq) = col_names
    prior_error <- counts_pnq %>%
        dplyr::group_by(.data$position, .data$nucleotide)%>%
        dplyr::summarise(counts=sum(.data$counts))%>%
        dplyr::ungroup() %>%
        dplyr::group_by(.data$position)%>%
        dplyr::mutate(prior.error = .data$counts/sum(.data$counts)) %>%
        dplyr::ungroup()%>%
        dplyr::select(.data$position, .data$nucleotide, .data$prior.error)

    if(save){utils::write.table(prior_error, file = output_file, row.names = FALSE)}
    return(prior_error)
}
