#' Train SINGLE model
#'
#' Main function to train a SINGLE model in a set of reads of a reference / wild type sequence. To get the input data you will need to run before a minimap2 alignment and samtools counts.
#'
#' @param file_ref_pnq File containing the counts per position returned by samtools mpileup
#' @param output_prefix String. Prefix for output files
#' @param refseq_fasta Fasta file containing reference sequence
#' @param rates.matrix Mutation rate matrix: 4x5 matrix, each row/col representing a nucleotide (col adds deletion), and the values is the mutational rate from row to col.
#' @param mean.n.mutations Mean number of mutations expected (one number).
#' @param pos_start Numeric. Position to start analyzing, counting starts from 1 and it refers to reference used for minimap2 alignment.
#' @param pos_end  Numeric. Position to stop analyzing, counting starts from 1 and it refers to reference used for minimap2 alignment.
#' @param verbose Logical.
#' @param save_partial Logical. Should partial results be saved in files?
#' @param save_final Logical. Should final fits be saved in a file?
#' @return Creates file output_prefix_single_results.txt with SINGLE training results.
#' @details Before running single_train_function you have to align your INPUT data to a REFERENCE using minimap2 and count the nucleotides per position using samtools using these lines:
#'
#'\code{minimap2 -ax map-ont --sam-hit-only  REFERENCE.fasta INPUT.fastq >ALIGNMENT.sam}
#'
#'\code{samtools view -S -b ALIGNMENT.sam > ALIGNMENT.bam}
#'
#'\code{samtools sort ALIGNMENT.bam -o ALIGNMENT.sorted.bam }
#'
#'\code{samtools mpileup -Q 0 ALIGNMENT.sorted.bam > COUNTS.txt}
#' @examples
#' pos_start=1
#' pos_end = 10
#' file_ref_pnq = system.file("extdata", "REF_READS_PNQ.txt", package = "single")
#' refseq_file = system.file("extdata", "ref_seq.fasta", package = "single")
#' output_prefix = tempfile()
#' fits_table <- single_train(file_ref_pnq,output_prefix,refseq_file,
#'                            mutation_rate,5,pos_start,pos_end, verbose = FALSE)
#' print(head(fits_table))
#' @export single_train
single_train    <- function(file_ref_pnq,output_prefix,refseq_fasta,
                                rates.matrix=NULL,mean.n.mutations=NULL,
                                pos_start=NULL,pos_end=NULL,
                                verbose=TRUE, save_partial=FALSE,
                                save_final=TRUE){
    options(dplyr.summarise.inform = FALSE)

    ### Verify inputs
    if(!is.character(file_ref_pnq)|length(file_ref_pnq)>1   ){
        stop("file_ref_pnq: must be a character of length 1")
    }
    if(!is.character(output_prefix)|length(output_prefix)>1       ){
        stop("output_prefix: must be a character of length 1")
    }
    if(!is.character(refseq_fasta)|length(refseq_fasta)>1){
        stop("refseq_fasta: must be a character of length 1")
    }

    ref_seq <- load_ref_seq(refseq_fasta)
    if(is.null(pos_start)     | !is.numeric(pos_start)){
        pos_start <- 1; warning("single_train: pos_start set to 1\n")
    }
    if(is.null(pos_end)       | !is.numeric(pos_end)  ){
        pos_end <- length(ref_seq);
        warning("single_train: pos_end set to length(ref_seq):", pos_end,"\n")
    }

    if(is.null(rates.matrix)){rates.matrix=mutation_rate; warning("single_train: rates matrix set to default mutation_rate\n")}
    if(is.null(colnames(rates.matrix))){stop('single_train: Missing colnames of rates.matrix')}
    if(is.null(rownames(rates.matrix))){stop('single_train: Missing rownames of rates.matrix')}
    if(is.null(mean.n.mutations)){stop('single_train: Missing mean.n.mutations')}

    ### Names of output and auxiliary files
    outfile_counts_nucleotides <- paste0(output_prefix,"_countsPNQ.txt")
    outfile_prior_errors       <- paste0(output_prefix,"_single_prob_prior_errors.txt")
    outfile_prior_mutations    <- paste0(output_prefix,"_single_prob_prior_mutations.txt")
    outfile_fits               <- paste0(output_prefix,"_single_fit.txt")
    outfile_data               <- paste0(output_prefix,"_single_data.txt")
    outfile_evaluation         <- paste0(output_prefix,"_single_results.txt")

    if(verbose){cat("single_train\n")}
    ### Parse counts of nucleotide per position and qscore
    if(verbose){cat("\t Parsing counts file \n")}
    counts_pnq <- parse_countspnq(input_file=file_ref_pnq,
                                    output_file=outfile_counts_nucleotides,
                                    pos_start=pos_start,pos_end=pos_end,
                                    save=save_partial)

    # Calculate priors errors
    if(verbose){cat("\t Computing p_prior-error \n")}
    priors_errors <- prior_errors(counts_pnq=counts_pnq,
                                    output_file=outfile_prior_errors,
                                    ref_seq=ref_seq,save=save_partial)

    ## Compute information of priors probabilities
    if(verbose){cat("\t Computing p_prior-right\n")}
    prior_mutation     <- prior_mutations(rates.matrix = rates.matrix,
                                            mean.n.mut = mean.n.mutations,
                                            ref_seq = ref_seq,
                                            save = save_partial,
                                            output_file=outfile_prior_mutations)

    #Fit wt
    if(verbose){cat("\t Fitting\n")}
    fits <- fit_logregr(counts_pnq = counts_pnq,
                        ref_seq=ref_seq,
                        prior_error = priors_errors,
                        prior_mutation = prior_mutation,
                        save=save_partial,
                        output_file_fits=outfile_fits,
                        output_file_data=outfile_data,
                        verbose = verbose)

    if(verbose){cat("\t Evaluating \n")}
    evaluated_fits <- evaluate_fits(pos_range = c(1,pos_end-pos_start+1),
                    q_range = c(1,50),
                    output_file = outfile_evaluation,
                    data_fits = fits,
                    ref_seq = ref_seq, verbose = verbose, save=save_final)

    return(evaluated_fits)
}
