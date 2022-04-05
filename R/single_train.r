#' Train SINGLE model
#'
#' Main function to train a SINGLE model in a set of reads of a reference / wild type sequence. To get the input data you will need to run before a minimap2 alignment and samtools counts.
#'
#' @param bamfile File containing the counts per position returned by samtools mpileup
#' @param output String. Prefix for output files
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
#' refseq_fasta<- system.file("extdata", "ref_seq.fasta", package = "single")
#' train_reads_example <- system.file("extdata", "train_seqs_500.sorted.bam",
#'                                    package = "single")
#' train <- single_train(bamfile=train_reads_example,
#'                    refseq_fasta=refseq_fasta,
#'                    rates.matrix=mutation_rate,mean.n.mutations=5.4,
#'                    pos_start=1,pos_end=10)
#' print(head(train))
#' @export single_train
#' @import Rsamtools
#' @importFrom Biostrings width readDNAStringSet
single_train    <- function(bamfile,
                            output="results",
                            refseq_fasta,
                            rates.matrix=NULL,mean.n.mutations=NULL,
                            pos_start=NULL,pos_end=NULL,
                            verbose=TRUE, save_partial=FALSE,
                            save_final=FALSE){

    options(dplyr.summarise.inform = FALSE)
    ref_seq <- Biostrings::readDNAStringSet(refseq_fasta)
    if(is.null(pos_start)|!is.numeric(pos_start)){
        pos_start <- 1;
        warning("single_train: pos_start set to 1\n")
    }
    if(is.null(pos_end)|!is.numeric(pos_end)){
        pos_end <- width(ref_seq);
        warning("single_train: pos_end set to length(ref_seq):", pos_end,"\n")
    }
    if(is.null(rates.matrix)){
        rates.matrix=mutation_rate;
        warning("single_train: rates matrix set to default mutation_rate\n")
    }
    if(is.null(colnames(rates.matrix))){
        stop('Missing colnames of rates.matrix')
    }
    if(is.null(rownames(rates.matrix))){
        stop('Missing rownames of rates.matrix')
    }
    if(is.null(mean.n.mutations)){
        stop('single_train: Missing mean.n.mutations')
    }

    ### Names of output and auxiliary files
    outfile_prior_errors    <- paste0(output,"_SINGLe_p_prior_errors.txt")
    outfile_prior_mutations <- paste0(output,"_SINGLe_p_prior_mutations.txt")
    outfile_fits            <- paste0(output,"_SINGLe_fits.txt")
    outfile_data            <- paste0(output,"_SINGLe_data.txt")
    outfile_evaluation      <- paste0(output,"_SINGLe_results.txt")

    if(verbose){message("single_train\n")}

    ### Parse counts of nucleotide per position and qscore
    if(verbose){message("\t pileup by quality score \n")}

    counts_pnq <- pileup_by_QUAL(bam_file=bamfile,
                    pos_start=pos_start,
                    pos_end=pos_end)

    # Calculate priors errors
    if(verbose){message("\t p_prior-error \n")}
    p_prior_errors <- p_prior_errors(counts_pnq=counts_pnq,
                        output_file=outfile_prior_errors,
                        save=save_partial)

    ## Compute information of priors probabilities
    if(verbose){message("\t p_prior-right\n")}
    p_prior_mutations   <- p_prior_mutations(rates.matrix = rates.matrix,
                            mean.n.mut = mean.n.mutations,
                            ref_seq = ref_seq,
                            save = save_partial,
                            output_file=outfile_prior_mutations)

    #Fit wt
    if(verbose){message("\t Fitting\n")}
    fits <- fit_logregr(counts_pnq = counts_pnq,
                        ref_seq=ref_seq,
                        p_prior_errors = p_prior_errors,
                        p_prior_mutations = p_prior_mutations,
                        save=save_partial,
                        output_file_fits=outfile_fits,
                        output_file_data=outfile_data,
                        verbose = verbose)

    if(verbose){message("\t Evaluating \n")}
    evaluated_fits <- evaluate_fits(pos_range = c(1,pos_end-pos_start+1),
                        q_range = c(1,50),
                        output_file = outfile_evaluation,
                        data_fits = fits,
                        ref_seq = ref_seq,
                        verbose = verbose, save=save_final)

    return(evaluated_fits)
}
