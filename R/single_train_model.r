#' Train SINGLE model
#'
#' Main function to train a SINGLE model in a set of reads of a reference / wild type sequence. To get the input data you will need to run before a minimap2 alignment and samtools counts.
#'
#' @param input_counts_file File containing the counts per position returned by samtools mpileup
#' @param output_prefix String. Prefix for output files
#' @param reference_sequence_fastafile Fasta file containing reference sequence
#' @param rates.matrix Mutation rate matrix: 4x5 matrix, each row/col representing a nucleotide (col adds deletion), and the values is the mutational rate from row to col.
#' @param mean.n.mutations Mean number of mutations expected (one number).
#' @param pos_start Numeric. Position to start analyzing, counting starts from 1 and it refers to reference used for minimap2 alignment.
#' @param pos_end  Numeric. Position to stop analyzing, counting starts from 1 and it refers to reference used for minimap2 alignment.
#' @param verbose Logical.
#' @param save_partial Logical. Should partial results be saved in files?
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
#' input_counts_file = system.file("extdata", "example_train_countsPNQ.txt", package = "single")
#' ref_seq_file = system.file("extdata", "ref_seq.fasta", package = "single")
#' outfile_ex = "example"
#' single_train_model(input_counts_file,outfile_ex,ref_seq_file, mutation_rate,5,5,10)
#' unlink(paste0(outfile_ex,"_single_results.txt"))
#' @export single_train_model
single_train_model    <- function(input_counts_file,output_prefix,reference_sequence_fastafile,
                                rates.matrix=NULL,mean.n.mutations=NULL,pos_start=NULL,pos_end=NULL,
                                verbose=TRUE, save_partial=FALSE){
    options(dplyr.summarise.inform = FALSE)

    ### Verify inputs
    if(!is.character(input_counts_file)      | length(input_counts_file)>1   ){stop("single_train_model: Wrong class input_counts_file: must be a character of length 1")}
    if(!is.character(output_prefix)          | length(output_prefix)>1       ){stop("single_train_model: Wrong class output_prefix: must be a character of length 1")}
    if(!is.character(reference_sequence_fastafile)| length(reference_sequence_fastafile)>1){stop("single_train_model: Wrong class reference_sequence_fastafile: must be a character of length 1")}

    reference_sequence <- load_reference_sequence(reference_sequence_fastafile)
    if(is.null(pos_start)     | !is.numeric(pos_start)){pos_start <- 1; warning("single_train_model: pos_start set to 1\n")}
    if(is.null(pos_end)       | !is.numeric(pos_end)  ){pos_end <- length(reference_sequence); warning("single_train_model; pos_end set to length(reference_sequence):", pos_end,"\n")}

    if(is.null(rates.matrix)){rates.matrix=mutation_rate; warning("single_train_model: rates matrix set to default mutation_rate\n")}
    if(is.null(colnames(rates.matrix))){stop('single_train_model: Missing colnames of rates.matrix')}
    if(is.null(rownames(rates.matrix))){stop('single_train_model: Missing rownames of rates.matrix')}
    if(is.null(mean.n.mutations)){stop('single_train_model: Missing mean.n.mutations')}

    ### Names of output and auxiliary files
    outfile_counts_nucleotides <- paste0(output_prefix,"_countsPNQ.txt")
    outfile_prior_errors       <- paste0(output_prefix,"_single_prob_prior_errors.txt")
    outfile_prior_mutations    <- paste0(output_prefix,"_single_prob_prior_mutations.txt")
    outfile_fits               <- paste0(output_prefix,"_single_fit.txt")
    outfile_data               <- paste0(output_prefix,"_single_data.txt")
    outfile_evaluation         <- paste0(output_prefix,"_single_results.txt")

    if(verbose){cat("single_train_model\n")}
    ### Parse counts of nucleotide per position and qscore
    if(verbose){cat("\t Parsing counts file \n")}
    counts_pnq <- parse_nucleotides_per_qscore(input_file=input_counts_file,output_file=outfile_counts_nucleotides,
                                                pos_start=pos_start,pos_end=pos_end, save=save_partial)

    # Calculate priors errors
    if(verbose){cat("\t Computing p_prior-error \n")}
    priors_errors <- calculate_prior_errors(counts_pnq=counts_pnq,output_file=outfile_prior_errors,
                                            reference_sequence=reference_sequence,save=save_partial)

    ## Compute information of priors probabilities
    if(verbose){cat("\t Computing p_prior-right\n")}
    prior_mutation     <- calculate_prior_mutations(rates.matrix = rates.matrix, mean.n.mut = mean.n.mutations,
                                                reference_sequence = reference_sequence,
                                                save = save_partial, output_file=outfile_prior_mutations)

    #Fit wt
    if(verbose){cat("\t Fitting\n")}
    fits <- fit_logistic_regression(counts_pnq = counts_pnq,
                                    reference_sequence=reference_sequence,
                                    prior_error = priors_errors,
                                    prior_mutation = prior_mutation,
                                    save=save_partial,
                                    output_file_fits=outfile_fits,
                                    output_file_data=outfile_data,
                                    verbose = verbose)

    if(verbose){cat("\t Evaluating \n")}
    evaluate_fits(pos_range = c(1,pos_end-pos_start+1),q_range = c(1,50),output_file = outfile_evaluation,data_fits = fits,
                reference_sequence = reference_sequence, verbose = verbose)

    return(0)
}
