#' Evaluate SINGLE model
#'
#' Main function to evaluate a gene library using a SINGLE model.
#'
#' @param input_sam_file File containing the counts per position returned by samtools mpileup
#' @param single_fits Results of the SINGLE model as returned by single_train(). It can be either the output data.frame or the saved file.
#' @param output_prefix String. Prefix for output files
#' @param refseq_fasta Fasta file containing reference sequence
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
#' input_sam_file = system.file("extdata", "LIB_READS.sam", package = "single")
#' single_fits = system.file("extdata", "single_fits_eval_ex.txt", package="single")
#' ref_seq_file = system.file("extdata", "ref_seq.fasta", package = "single")
#' outfile = tempfile("example")
#' single_evaluate(input_sam_file,single_fits,outfile,ref_seq_file,
#'     pos_start=1,pos_end=100,verbose=TRUE)
#' @export single_evaluate
single_evaluate <- function(input_sam_file,
                            single_fits,
                            output_prefix,
                            refseq_fasta,
                            pos_start=NULL,pos_end=NULL,
                            verbose=TRUE,
                            gaps_weights="none"){

    ## Verify inputs
    if(!is.character(input_sam_file)| length(input_sam_file)>1){ stop("input_sam_file: must be a character of length 1")}
    if( !(is.character(single_fits) & length(single_fits)==1) & !(is.data.frame(single_fits)&& ncol(single_fits)==5)){
        stop("single_fits: must be a character of length 1 ir a data.frame")
    }
    if(!is.character(output_prefix) | length(output_prefix)>1){  stop("single_evaluate: Wrong class output_prefix: must be a character of length 1")}
    if(!is.character(refseq_fasta)| length(refseq_fasta)>1){stop("single_evaluate: Wrong class refseq_fasta: must be a character of length 1")}
    ref_seq <- load_ref_seq(refseq_fasta)
    if(is.null(pos_start) | !is.numeric(pos_start)){pos_start <- 1; warning("pos_start set to 1\n")}
    if(is.null(pos_end) | !is.numeric(pos_end)){pos_end <- length(ref_seq); warning("pos_end set to length(ref_seq)\n")}
    if(length(gaps_weights)>1 | !gaps_weights %in%c("minimum","none","mean")){stop("gaps_weight should be minimum, none, or mean")}

    ## Define outputs names
    output_txt_file   <- paste0(output_prefix,"_corrected.txt")

    #Load fitted values and priors
    if(verbose){cat("\nLoad fits data \n")}

    if(verbose){cat("\nEvaluate sequences and save results \n")}
    save_txt_file(input_sam_file=input_sam_file,
                    single_fits=single_fits,
                    output_txt_file=output_txt_file,
                    gaps_weights=gaps_weights,
                    ref_seq=ref_seq,
                    pos_start=pos_start, pos_end=pos_end,
                    comment.char = "@")

    return(0)
}
