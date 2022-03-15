#' Evaluate SINGLE fits
#'
#' Evaluates single fits for all positions, nucleotides and Qscores in the given ranges.
#' @param pos_range Numeric vector. Positions to evaluate.
#' @param q_range Numeric vector. Qscores to evaluate.
#' @param output_file File name for output, if save=TRUE.
#' @param data_fits Data.frame with columns position nucleotide slope intercept as the one returned by fit_logregr
#' @param ref_seq Reference sequence: vector of characters, as returned by load_ref_seq
#' @param verbose Logical.
#' @param save Logical. Should results be saved in output_file?
#' @return data.frame with SINGLE fits evaluated for pos_range and q_range. Columns are: position, nucleotide, quality, p_SINGLe, isWT.
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom utils txtProgressBar setTxtProgressBar write.table
#' @export evaluate_fits
#' @examples
#' pos_range = seq_len(100)
#' q_range = seq(1,50)
#' ref_seq_file = system.file("extdata", "ref_seq.fasta", package = "single")
#' ref_seq = readDNAStringSet(ref_seq_file)
#' evaluated_fits <- evaluate_fits(pos_range = c(1,5),q_range = c(0,10),
#'                      data_fits = fits,ref_seq = ref_seqE)
evaluate_fits <- function(pos_range,q_range,output_file, data_fits, ref_seq,
                          verbose=FALSE, save=FALSE){
    t0 <- proc.time()
    ref_seq_char = strsplit(as.character(ref_seq),"")[[1]]

    frequencies <- expand.grid(seq(pos_range[1],pos_range[2]),bases,seq(q_range[1],q_range[2]), c("+","-")) %>%
        dplyr::rename(pos = .data$Var1,nucleotide=.data$Var2,QUAL=.data$Var3, strand=.data$Var4) %>%
        dplyr::arrange(.data$pos,.data$nucleotide,.data$QUAL) %>%
        dplyr::mutate(wt.base=ref_seq_char[.data$pos])     %>%          # Column with reference nucleotide
        dplyr::mutate(isWT = .data$nucleotide==.data$wt.base)         %>%          # Logical column indicating if read is reference or different nucleotide
        dplyr::mutate(p_SINGLe=NA)

    frequencies <- dplyr::left_join(frequencies,data_fits,by=c("pos","nucleotide","strand"))
    if(verbose){pb = utils::txtProgressBar(min = 0, max = nrow(frequencies), style = 3)}
    for(i in seq_len(nrow(frequencies))){
        if(verbose){utils::setTxtProgressBar(pb,i)}
        # if it is a wild type, skip
        if(frequencies$isWT[i]==1){ next() }
        #fit models
        frequencies$p_SINGLe [i] <- glm.predict.(x=frequencies$QUAL[i],
                                                             slope = frequencies$prior_slope[i],
                                                             intercept = frequencies$prior_intercept[i])
    }
    t1 <- proc.time()
    if(verbose){print(t1-t0)}

    frequencies$p_SINGLe[which(frequencies$isWT==TRUE)] <- 1-10^(-frequencies$QUAL[which(frequencies$isWT==TRUE)]/10)     # For wildtype keep original Qscore (or values that did not occurred)
    frequencies <- frequencies %>%
        dplyr::select(.data$pos,.data$strand, .data$nucleotide, .data$QUAL, .data$p_SINGLe, .data$isWT)
    #Save results
    if(save){
        utils::write.table(frequencies, file = output_file, row.names = FALSE)
    }
    return(frequencies)
}


