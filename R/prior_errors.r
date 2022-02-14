#' Computes prior probability of errors
#'
#' This is an auxiliary function in single package. It takes a data frame with counts by position, nucleotide and Qscore and it summarises it into proportion of nucleotide counts by position.
#' @param counts_pnq Data frame with columns position nucleoide quality counts, as returned by parse_countspnq
#' @param output_file File name for output, if save=TRUE.
#' @param ref_seq Reference sequence: vector of characters, as returned by load_ref_seq
#' @param save Logical. Should data be saved in a output_file?
#' @return Data frame with columns position nucleotide prior.error.
#' @importFrom utils write.table
#' @importFrom rlang .data
#' @export prior_errors
#' @examples
#' refseq_file = system.file("extdata", "ref_seq.fasta", package = "single")
#' ref_seq = load_ref_seq(refseq_file)
#' file_ref_pnq = system.file("extdata", "REF_READS_PNQ.txt", package = "single")
#' counts_pnq = parse_countspnq(input_file=file_ref_pnq,output_file=NA)
#' prior_errors(counts_pnq,ref_seq=ref_seq)
prior_errors         <- function(counts_pnq,output_file=NULL, ref_seq,save=FALSE){
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
