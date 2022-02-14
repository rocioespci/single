#' Evaluate SINGLE fits
#'
#' Evaluates single fits for all positions, nucleotides and Qscores in the given ranges.
#' @param pos_range Numeric vector. Positions to evaluate.
#' @param q_range Numeric vector. Qscores to evaluate.
#' @param output_file File name for output, if save=TRUE.
#' @param data_fits Data.frame with columns position nucleotide slope intercept as the one returned by fit_logistic_regression.
#' @param reference_sequence Reference sequence: vector of characters, as returned by load_reference_sequence
#' @param verbose Logical.
#' @return data.frame with SINGLE fits evaluated for pos_range and q_range. Columns are: position, nucleotide, quality, p_right_priors_model, isWT.
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom utils txtProgressBar setTxtProgressBar write.table
#' @export evaluate_fits
#' @examples
#' pos_range = seq_len(10)
#' q_range = seq(10,30)
#' ref_seq_file = system.file("extdata", "ref_seq.fasta", package = "single")
#' ref_seq = load_reference_sequence(ref_seq_file)
#' counts_pnq_file = system.file("extdata", "example_train_countsPNQ_sub.txt", package = "single")
#' counts_pnq = read.table(counts_pnq_file, header=TRUE)
#' prior_error = calculate_prior_errors(counts_pnq,reference_sequence=ref_seq)
#' prior_mutation = calculate_prior_mutations(mutation_rate, 3, ref_seq)
#' fits = fit_logistic_regression(counts_pnq, ref_seq,prior_error,prior_mutation)
#' outfile_ex = tempfile("evaluate_fits_example.txt")
#' evaluate_fits(pos_range,q_range,outfile_ex, fits,ref_seq, verbose=FALSE)
evaluate_fits                  <- function(pos_range,q_range,output_file, data_fits, reference_sequence, verbose=TRUE){
    t0 <- proc.time()
    frequencies <- expand.grid(seq(pos_range[1],pos_range[2]),bases,seq(q_range[1],q_range[2])) %>%
        dplyr::rename(position = .data$Var1,nucleotide=.data$Var2,quality=.data$Var3) %>%
        dplyr::arrange(.data$position,.data$nucleotide,.data$quality) %>%
        dplyr::mutate(wt.base=reference_sequence[.data$position])     %>%          # Column with reference nucleotide
        dplyr::mutate(isWT = .data$nucleotide==.data$wt.base)         %>%          # Logical column indicating if read is reference or different nucleotide
        dplyr::mutate(p_right_priors_model=NA)

    frequencies <- dplyr::left_join(frequencies,data_fits,by=c("position","nucleotide"))
    if(verbose){pb = utils::txtProgressBar(min = 0, max = nrow(frequencies), style = 3)}
    for(i in seq_len(nrow(frequencies))){
        if(verbose){utils::setTxtProgressBar(pb,i)}
        # if it is a wild type, skip
        if(frequencies$isWT[i]==1){ next() }
        #fit models
        frequencies$p_right_priors_model [i] <- glm.predict.(x=frequencies$quality[i],slope = frequencies$prior_slope[i], intercept = frequencies$prior_intercept[i])
    }
    t1 <- proc.time()
    if(verbose){print(t1-t0)}

    frequencies$p_right_priors_model[which(frequencies$isWT==TRUE)] <- 1-10^(-frequencies$quality[which(frequencies$isWT==TRUE)]/10)     # For wildtype keep original Qscore (or values that did not occurred)
    frequencies <- frequencies %>%
        dplyr::select(.data$position, .data$nucleotide, .data$quality, .data$p_right_priors_model, .data$isWT)
    #Save results
    utils::write.table(frequencies, file = output_file, row.names = FALSE)
    return(frequencies)
}


