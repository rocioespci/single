#' Fit SINGLE's logistic regression
#'
#' This is an auxiliary function in single package. It takes counts_pnq and for each position and nucleotide it fits SINGLE's logistic regression.
#' @param counts_pnq Data frame with columns position nucleoide quality counts, as returned by pileup_by_QUAL
#' @param ref_seq Reference sequence: vector of characters, as returned by load_ref_seq
#' @param p_prior_errors Data frame with columns position nucleotide prior.error, as the one returned by p_prior_errors().
#' @param p_prior_mutations Data frame with columns wt.base, nucleotide and p_mutation (probaility of mutation), as the one returned by p_prior_mutations().
#' @param save Logical. Should data be saved in a output_file?
#' @param output_file_fits File into which save the single fits if save=TRUE
#' @param output_file_data File into which save the fitted data if save=TRUE
#' @param verbose Logical.
#' @return data.frame with columns position, nucleotide, slope and intercept (of the sigmoidal regression).
#' @importFrom tidyr replace_na
#' @importFrom stats glm coefficients
#' @importFrom utils write.table
#' @importFrom rlang .data
#' @import dplyr
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export fit_logregr
#' @examples
#' pos_start <- 1
#' pos_end <- 10
#' refseq_fasta <- system.file("extdata", "ref_seq.fasta", package = "single")
#' train_reads_example <- system.file("extdata", "train_seqs_500.sorted.bam",
#'                                    package = "single")
#' counts_pnq <- pileup_by_QUAL(bam_file=train_reads_example,
#'     pos_start=pos_start,pos_end=pos_end)
#' p_prior_mutations <- p_prior_mutations(rates.matrix = mutation_rate,
#'     mean.n.mut = 5,ref_seq = ref_seq,save = F)
#' p_prior_errors <- p_prior_errors(counts_pnq=counts_pnq,save=FALSE)
#' fits <- fit_logregr(counts_pnq = counts_pnq,ref_seq=ref_seq,
#'     p_prior_errors = p_prior_errors,p_prior_mutations = p_prior_mutations,
#'     save=F)
#' head(fits)
fit_logregr <- function(counts_pnq,ref_seq,p_prior_errors,p_prior_mutations,
                save=FALSE, output_file_fits,output_file_data,verbose=FALSE){

    #Pre-editing data
    ref_seq_char = strsplit(as.character(ref_seq),"")[[1]]
    data <- dplyr::as_tibble(counts_pnq)%>%
        dplyr::mutate(wt.base = ref_seq_char[.data$pos])          # Add wildtype base

    ## Wildtype matrix: keep rows with wildtype reads
    data_wt <- data %>%
        dplyr::filter(.data$nucleotide==.data$wt.base)%>%
        dplyr::select(-.data$nucleotide) %>%
        dplyr::rename(counts.wt=count)

    ## Data with mutations (errors)
    data_mut <- data %>%                                             #start from data
        dplyr::filter(.data$nucleotide!=.data$wt.base)              #keep only mutations

    missing_position <- setdiff(seq_len(max(data$pos)) ,data_mut$pos)
    if(length(missing_position)>0){
        data.aux <- data.frame(pos=missing_position, nucleotide=NA, QUAL=20,count=NA,wt.base=NA) %>%
            mutate(wt.base=ref_seq_char[.data$pos]) %>%
            mutate(nucleotide=if_else(.data$wt.base=="A","C","A")) %>%
            mutate(strand="+")
        data_mut <- data_mut %>% rbind(data.aux)
        rm(data.aux)
    }

    data_mut_expansion <- data_mut %>%
        tidyr::expand(.data$pos,.data$nucleotide,.data$QUAL,.data$strand)
    data_mut <- data_mut %>%
        dplyr::full_join(data_mut_expansion,
                        by = c("pos", "nucleotide", "QUAL","strand"))%>%        # complete all combinations position - nucleotide - quality
        dplyr::mutate(wt.base = ref_seq_char[.data$pos]) %>%                    # fill wildtype base for missing values (new rows)
        dplyr::full_join(data_wt,   by=c("pos", "QUAL","wt.base","strand"))%>%   # add wildtype info in new columns
        dplyr::left_join(p_prior_mutations, by=c("wt.base", "nucleotide"))%>%      # add prior of being mutated
        dplyr::left_join(p_prior_errors,    by=c("pos", "nucleotide","strand"))%>%     # add prior of being an error
        dplyr::filter(.data$nucleotide!=.data$wt.base)%>%                                   # remove wildtype rows
        dplyr::mutate(count=tidyr::replace_na(.data$count,0))%>%                           # fill NA with 0
        dplyr::mutate(counts.wt=tidyr::replace_na(.data$counts.wt,0))%>%                     # fill NA with 0
        dplyr::mutate(p_prior_error=tidyr::replace_na(.data$p_prior_error,0))                    # fill NA with 0
    rm(data_mut_expansion)

    ## Count number of errors/good reads by position and nucleotide
    total_counts <- data_mut %>%
        dplyr::group_by(.data$pos, .data$nucleotide,.data$strand)%>%
        dplyr::summarise(total.counts.mut=sum(.data$count,na.rm=TRUE),
                        total.counts.wt=sum(.data$counts.wt,na.rm=TRUE))%>%
        dplyr::mutate(total.counts=.data$total.counts.mut+.data$total.counts.wt)
    data_mut <- dplyr::full_join(data_mut, total_counts, by=c("pos", "nucleotide","strand"))

    ## reweight counts by prior probabilities
    data_mut <- data_mut %>%
        dplyr::mutate(pc = .data$p_mutation / (.data$p_mutation+.data$p_prior_error),                            #prob of being correct
                        pi = .data$p_prior_error / (.data$p_mutation+.data$p_prior_error))%>%                        #prob of being an error
        dplyr::mutate(counts.scaled    = (.data$count / .data$total.counts.mut * .data$pi * .data$total.counts ),   #counts errors re-weighted
                    counts.wt.scaled = (.data$counts.wt / .data$total.counts.wt * .data$pc * .data$total.counts )) #counts wildtype re-weighted

    ## Fit data:
    data_fits <- data_mut %>%
        dplyr::select(.data$strand,.data$pos, .data$nucleotide)%>%
        dplyr::distinct(.data$strand,.data$pos, .data$nucleotide)%>%
        dplyr::arrange(.data$strand,.data$pos, .data$nucleotide)%>%
        dplyr::mutate(prior_slope=NA, prior_intercept=NA)
    n_i <- nrow(data_fits)
    if(verbose){cat("\n Fitting \n")}
    if(verbose){ pb = utils::txtProgressBar(min = 0, max = nrow(data_fits), style = 3) }
    for (i in seq_len(nrow(data_fits))){
        if(verbose){utils::setTxtProgressBar(pb,i)}
        #Keep data from this position & nucleotide
        aux_df <- data_mut %>%
            dplyr::filter(.data$pos==data_fits$pos[i] & .data$nucleotide ==data_fits$nucleotide[i] & .data$strand==data_fits$strand[i])%>%
            dplyr::select(.data$strand,.data$QUAL,.data$count, .data$counts.wt,.data$counts.scaled,.data$counts.wt.scaled) %>%
            dplyr::mutate(tot.scaled=.data$counts.scaled+.data$counts.wt.scaled) %>%
            dplyr::mutate(proportion.wt.scaled=.data$counts.wt.scaled/.data$tot.scaled)

        if(sum(aux_df$count)==0){
            data_fits$prior_slope[i]       <- NA
            data_fits$prior_intercept[i]   <- NA

            if(verbose){warning('Position ', data_fits$pos[i], data_fits$nucleotide[i], " has no data to be fitted.\n")}
            next()
        }
        qval  <- aux_df$QUAL
        yvals <- aux_df$proportion.wt.scaled

        #Corrected fit by prior data
        if(! (all(is.na(aux_df$counts.wt.scaled)) | all(is.na(aux_df$counts.scaled))) ){
            aux_prior_fit                <- stats::glm(yvals[!is.na(yvals)]~ qval[!is.na(yvals)],family = "quasibinomial")
            aux_prior_coefficients       <- stats::coefficients(aux_prior_fit)
            data_fits$prior_slope[i]     <- aux_prior_coefficients[2]
            data_fits$prior_intercept[i] <- aux_prior_coefficients[1]
        }
    }
    # SAVE RESULTS
    if(save){
        utils::write.table(data_fits, file=output_file_fits, quote = FALSE, row.names = FALSE)
        utils::write.table(data_mut,file=output_file_data, quote = FALSE, row.names = FALSE)
    }
    return(data_fits)
}
