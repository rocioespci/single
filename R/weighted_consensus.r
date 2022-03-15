#' Compute consensus sequence
#'
#' This is an auxiliary function in single package. It computes consensus from a data.frame as the one returned by single_evaluate()
#' @param df data.frame with the columns: nucleotide, probability, position
#' @param cutoff_prob Numeric. Nucleotides with probability below this number will be removed from consensus computation.
#' @return Character vector, consensus sequence
#' @importFrom rlang .data
#' @importFrom dplyr %>% select rename
#' @importFrom Biostrings DNAString
#' @export weighted_consensus
#' @examples
#' data_barcode = data.frame(
#'  nucleotide = unlist(sapply(as.character(corrected_reads),strsplit, split="")),
#'  p_SINGLe=unlist(1-as(quality(corrected_reads),"NumericList")),
#'  pos=rep(1:width(corrected_reads[1]),length(corrected_reads)))
#' weighted_consensus(df = data_barcode, cutoff_prob = 0.9)
#' weighted_consensus(df = data_barcode, cutoff_prob = 0.999)
weighted_consensus       <- function(df,cutoff_prob=0.2){
    if(cutoff_prob>1 | cutoff_prob<0){stop('cutoff must be between 0 and 1')}

    checksequence <- setdiff(unique(df[,1]),bases)
    if(length(checksequence)>0){
        stop('consensus_seq: df first column should be only the nucleotides')
    }
    checkprob <- range(df[,2], na.rm=TRUE)
    if(!is.numeric(checkprob) || any(checkprob<0,na.rm = TRUE)|| any(checkprob>1,na.rm=TRUE)){
        stop('consensus_seq: df second column should be probabilities')
    }
    if(!is.numeric(df[,3]) & !is.integer(df[,3]) ){
        stop('consensus_seq: third column should be position')
    }
    colnames(df) <- c("nucleotide","probability","pos")
    df.counts <- df %>%
        dplyr::filter(.data$probability>cutoff_prob)     %>%
        dplyr::group_by(.data$pos, .data$nucleotide)     %>%
        dplyr::summarise(counts=sum(.data$probability))  %>%
        dplyr::ungroup()                                 %>%
        dplyr::group_by(.data$pos)                       %>%
        dplyr::filter(.data$counts==max(.data$counts))

    #If there are 2 bases which are max, I keep the first one (some bias towards A/C..)
    if(any(duplicated(df.counts$pos))){
        df.counts <- df.counts[!duplicated(df.counts$pos),]
    }

    #Complete all positions
    missing.positions <- setdiff(seq_len(max(df$pos)), df.counts$pos)
    if(length(missing.positions)>0){
        aux_df <- data.frame(pos=missing.positions,
                    nucleotide=factor("N", levels=levels(df.counts$nucleotide)),
                    counts=0)
        df.counts <- df.counts %>%
            dplyr::bind_rows(aux_df) %>%
            dplyr::arrange(.data$pos)
    }

    # return(as.character(df.counts$nucleotide))
    return(Biostrings::DNAString(paste0(df.counts$nucleotide, collapse = "")))
}
