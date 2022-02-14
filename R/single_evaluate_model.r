#' Evaluate SINGLE model
#'
#' Main function to evaluate a gene library using a SINGLE model.
#'
#' @param input_sam_file File containing the counts per position returned by samtools mpileup
#' @param input_fits_file File containing the results of the SINGLE model. It will be saved when you run single_train_model.
#' @param output_prefix String. Prefix for output files
#' @param reference_sequence_fastafile Fasta file containing reference sequence
#' @param pos_start Numeric. Position to start analyzing, counting starts from 1 and it refers to reference used for minimap2 alignment.
#' @param pos_end  Numeric. Position to stop analyzing, counting starts from 1 and it refers to reference used for minimap2 alignment.
#' @param gaps_weights One of "minimum","none","mean". How to assign qscores to deletions.
#' @param verbose Logical
#' @return Creates file output_prefix_corrected.txt with the Qscores re-scaled by SINGLE. Columns are SeqID position nucleotide isWT original_quality p_right_priors_model.
#' @details Before running single_evaluate_function you have to align your INPUT data to a REFERENCE using minimap2 and count the nucleotides per position using samtools using these lines:
#'
#'\code{minimap2 -ax map-ont --sam-hit-only  REFERENCE.fasta INPUT.fastq >ALIGNMENT.sam}
#'
#'\code{samtools view -S -b ALIGNMENT.sam > ALIGNMENT.bam}
#'
#'\code{samtools sort ALIGNMENT.bam -o ALIGNMENT.sorted.bam }
#'
#'\code{samtools mpileup -Q 0 ALIGNMENT.sorted.bam > COUNTS.txt}
#' @examples
#' input_sam_file = system.file("extdata", "example_sequences.sam", package = "single")
#' input_fits_file = system.file("extdata", "example_single_results.txt", package="single")
#' ref_seq_file = system.file("extdata", "ref_seq.fasta", package = "single")
#' outfile_ex = tempfile("example")
#' single_evaluate_model(input_sam_file,input_fits_file,outfile_ex,ref_seq_file, verbose=TRUE)
#' unlink(paste0(outfile_ex,"_single_results.txt"))
#' @export single_evaluate_model
single_evaluate_model <- function(input_sam_file,
                                  input_fits_file,
                                  output_prefix,
                                  reference_sequence_fastafile,
                                  pos_start=NULL,pos_end=NULL,
                                  verbose=TRUE,
                                  gaps_weights="none"){

    ## Verify inputs
    if(!is.character(input_sam_file)| length(input_sam_file)>1){ stop("single_evaluate_model: Wrong class input_sam_file: must be a character of length 1")}
    if(!is.character(input_fits_file)| length(input_fits_file)>1){stop("single_evaluate_model: Wrong class input_fits_file: must be a character of length 1")}
    if(!is.character(output_prefix) | length(output_prefix)>1){  stop("single_evaluate_model: Wrong class output_prefix: must be a character of length 1")}
    if(!is.character(reference_sequence_fastafile)| length(reference_sequence_fastafile)>1){stop("single_evaluate_model: Wrong class reference_sequence_fastafile: must be a character of length 1")}
    reference_sequence <- load_reference_sequence(reference_sequence_fastafile)
    if(is.null(pos_start) | !is.numeric(pos_start)){pos_start <- 1; warning("pos_start set to 1\n")}
    if(is.null(pos_end) | !is.numeric(pos_end)){pos_end <- length(reference_sequence); warning("pos_end set to length(reference_sequence)\n")}
    if(length(gaps_weights)>1 | !gaps_weights %in%c("minimum","none","mean")){stop("gaps_weight should be minimum, none, or mean")}

    ## Define outputs names
    output_txt_file   <- paste0(output_prefix,"_corrected.txt")

    #Load fitted values and priors
    if(verbose){cat("\nLoad fits data \n")}

    if(verbose){cat("\nEvaluate sequences and save results \n")}
    save_txt_file(input_sam_file=input_sam_file,
                  aux_fitted_file=input_fits_file,
                  output_txt_file=output_txt_file,
                  gaps_weights=gaps_weights,
                  reference_sequence=reference_sequence,
                  pos_start=pos_start, pos_end=pos_end,
                  comment.char = "@")

  return(0)
}
