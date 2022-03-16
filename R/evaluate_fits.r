#' Evaluate SINGLE fits
#'
#' Evaluates SINGLe for pos, nucleotides and QUAL in the given ranges.
#' @param pos_range Numeric vector. Positions to evaluate.
#' @param q_range Numeric vector. QUAL to evaluate.
#' @param output_file File name for output, if save=TRUE.
#' @param data_fits Data.frame with columns position nucleotide slope intercept as the one returned by fit_logregr
#' @param ref_seq DNAStringSet containing the true reference sequence.
#' @param verbose Logical.
#' @param save Logical. Should results be saved in output_file?
#' @return data.frame with SINGLE fits evaluated for pos_range and q_range.
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom utils txtProgressBar setTxtProgressBar write.table
#' @importFrom Biostrings readDNAStringSet
#' @export evaluate_fits
#' @examples
#' pos_range = seq_len(100)
#' q_range = seq(1,50)
#' ref_seq_file = system.file("extdata", "ref_seq.fasta", package = "single")
#' ref_seq = Biostrings::readDNAStringSet(ref_seq_file)
#' evaluated_fits <- evaluate_fits(pos_range = c(1,5),q_range = c(0,10),
#'                      data_fits = fits,ref_seq = ref_seqE)
evaluate_fits <- function(pos_range,q_range,output_file, data_fits, ref_seq,
                            verbose=FALSE, save=FALSE){
    t0 <- proc.time()
    ref_seq_char = strsplit(as.character(ref_seq),"")[[1]]

    df <- expand.grid(
            seq(pos_range[1],pos_range[2]),
            bases,seq(q_range[1],q_range[2]), c("+","-")) %>%
        dplyr::rename(pos = .data$Var1,nucleotide=.data$Var2,
                            QUAL=.data$Var3, strand=.data$Var4) %>%
        dplyr::arrange(.data$pos,.data$nucleotide,.data$QUAL) %>%
        dplyr::mutate(wt.base=ref_seq_char[.data$pos])%>%
        dplyr::mutate(isWT = .data$nucleotide==.data$wt.base) %>%
        dplyr::mutate(p_SINGLe=NA)

    df <- dplyr::left_join(df,data_fits,
                                    by=c("pos","nucleotide","strand"))
    if(verbose){p =utils::txtProgressBar(min =0,max=nrow(df),style = 3)}
    for(i in seq_len(nrow(df))){
        if(verbose){utils::setTxtProgressBar(pb,i)}
        # if it is a wild type, skip
        if(df$isWT[i]==1){ next() }
        #fit models
        df$p_SINGLe [i] <- glm.predict.(x=df$QUAL[i],
                                    slope = df$prior_slope[i],
                                    intercept = df$prior_intercept[i])
    }
    t1 <- proc.time()
    if(verbose){print(t1-t0)}
    # For wildtype keep original Qscore (or values that did not occurred):
    df$p_SINGLe[which(df$isWT)] <- 1-10^(-df$QUAL[which(df$isWT)]/10)
    df <- df %>%
        dplyr::select(.data$pos,.data$strand, .data$nucleotide,
                        .data$QUAL, .data$p_SINGLe, .data$isWT)
    #Save results
    if(save){
        utils::write.table(df, file = output_file, row.names = FALSE)
    }
    return(df)
}


