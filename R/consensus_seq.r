#' Compute consensus sequence
#'
#' This is an auxiliary function in single package. It computes consensus from a data.frame as the one returned by single_evaluate()
#' @param df data.frame with the columns: nucleotide, probability, position
#' @param cutoff_prob Numeric. Nucleotides with probability below this number will be removed from consensus computation.
#' @return Character vector, consensus sequence
#' @importFrom rlang .data
#' @importFrom dplyr %>% select rename
#' @export consensus_seq
#' @examples
#' require(dplyr)
#' seqs = system.file("extdata", "single_reads_ex_corrected.txt", package="single")
#' seqs_table= read.table(seqs, header=TRUE,nrows=100*10+1)
#' seqs_table = seqs_table %>%
#'     select(nucleotide,p_right_priors_model,position) %>%
#'     rename(probability=p_right_priors_model)
#' consensus_seq(seqs_table)
consensus_seq       <- function(df,cutoff_prob=0.2){
    if(cutoff_prob>1 | cutoff_prob<0){stop('cutoff must be between 0 and 1')}

    checksequence <- setdiff(bases,unique(df[,1]))
    if(length(checksequence)>0){ stop('consensus sequence v2: df first column should be only the nucleotides')}else{
        colnames(df)[1] <- "nucleotide"
    }
    checkprobabilities <- range(df[,2], na.rm=TRUE)
    if(!is.numeric(checkprobabilities) || any(checkprobabilities<0 , na.rm = TRUE)|| any(checkprobabilities>1,na.rm=TRUE)){
        stop('consensus sequence v2: df second column should be probabilities')}else{
            colnames(df)[2] <- "probability"
        }
    if(!is.numeric(df[,3]) & !is.integer(df[,3]) ){
        stop('consensus sequence v2: third column should be position')
    }else{
        colnames(df)[3] <- "position"
    }

    df.counts <- df %>%
        dplyr::filter(.data$probability>cutoff_prob)     %>%
        dplyr::group_by(.data$position, .data$nucleotide)      %>%
        dplyr::summarise(counts=sum(.data$probability))  %>%
        dplyr::ungroup()                           %>%
        dplyr::group_by(.data$position)                  %>%
        dplyr::filter(.data$counts==max(.data$counts))

    #If there are 2 bases which are max, I keep the first one (some bias towards A/C..)
    if(any(duplicated(df.counts$position))){  df.counts <- df.counts[!duplicated(df.counts$position),]  }

    #Complete all positions
    missing.positions <- setdiff(seq_len(max(df$position)), df.counts$position)
    if(length(missing.positions)>0){
        aux_df <- data.frame(position=missing.positions,
                            nucleotide=factor("N", levels=levels(df.counts$nucleotide)),
                            counts=0)
        df.counts <- df.counts %>%
            dplyr::bind_rows(aux_df) %>%
            dplyr::arrange(.data$position)
    }

    # return(as.character(df.counts$nucleotide))
    return(Biostrings::DNAString(paste0(df.counts$nucleotide, collapse = "")))
}
