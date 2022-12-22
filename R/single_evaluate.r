#' Evaluate SINGLE model
#'
#' Main function to evaluate a gene library using a SINGLE model.
#'
#' @param bamfile File containing the counts per position returned by samtools mpileup
#' @param single_fits Results of the SINGLE model as returned by single_train(). It can be either the output data.frame or the saved file.
#' @param refseq_fasta Fasta file containing reference sequence
#' @param pos_start Numeric. Position to start analyzing, counting starts from 1 and it refers to reference used for minimap2 alignment.
#' @param pos_end  Numeric. Position to stop analyzing, counting starts from 1 and it refers to reference used for minimap2 alignment.
#' @param gaps_weights One of "minimum","none","mean". How to assign qscores to deletions.
#' @param save Logical. Should data be saved in a output_file?
#' @param output_file File name for output, if save=TRUE.
#' @param verbose Logical
#' @param save_original_scores Logical. Should original Qscores be saved? If TRUE, they are stored in file whose name finishes in _original.fastq
#' @return Creates file output_prefix_corrected.txt with the Qscores re-scaled by SINGLE. Columns are SeqID position nucleotide isWT original_quality p_SINGLe
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
#' refseq_fasta = system.file("extdata", "ref_seq_10bases.fasta", package = "single")
#' train_file <- system.file("extdata", "train_example.txt", package = "single")
#' train <- read.table(train_file, header=TRUE)
#' test_reads_example <- system.file("extdata", "test_sequences.sorted.bam",
#'    package = "single")
#' corrected_reads <- single_evaluate(bamfile = test_reads_example,
#'                  single_fits = train,refseq_fasta = refseq_fasta,
#'                  pos_start=1,pos_end=10,gaps_weights = "minimum")
#' corrected_reads
#' @export single_evaluate
#' @importFrom Biostrings QualityScaledDNAStringSet writeQualityScaledXStringSet PhredQuality readDNAStringSet vmatchPattern replaceAt extractAt compareStrings width
#' @importFrom dplyr %>% left_join select
#' @importFrom IRanges IRanges
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stringr str_locate_all
#' @importFrom Rsamtools BamFile scanBam
#' @importFrom GenomicAlignments sequenceLayer
#' @importFrom Biostrings subseq
#' @importFrom BiocGenerics start
#' @importFrom methods as
single_evaluate <- function(bamfile, single_fits,
                            refseq_fasta,
                            pos_start=NULL, pos_end=NULL,
                            gaps_weights,
                            save=FALSE,output_file,
                            verbose=FALSE,
                            save_original_scores=FALSE){

    #Verify inputs
    if(is.null(pos_start)){pos_start=1;warning("save_txt_file: pos_start set to 1")}
    if(is.null(pos_end)){pos_end=width(ref_seq);
        warning("save_txt_file: pos_end set to ",length(ref_seq)+pos_start,"\n")}
    if(!file.exists(bamfile)){stop("save_txt_file: file ", bamfile," does not exist")}

    ref_seq <- Biostrings::readDNAStringSet(refseq_fasta)

    # Load data and define general parameters
    if(is.character(single_fits)){
        if(!file.exists(single_fits)){
            stop("save_txt_file: file ", single_fits," does not exist")
        }
        qtable = utils::read.table(single_fits, header=TRUE)
    }else{
        qtable = single_fits
    }

    bf <- Rsamtools::BamFile(bamfile)
    reads <- Rsamtools::scanBam(bf)
    #keep sequences that start at least at pos_start
    index <- which(reads[[1]]$pos<=pos_start)
    reads[[1]] <- sapply(reads[[1]], function(x){x[index]})

    reads_aligned <- GenomicAlignments::sequenceLayer(reads[[1]]$seq,
                                reads[[1]]$cigar,to = "reference")
    names(reads_aligned) <- reads[[1]]$qname
    scores_aligned <- GenomicAlignments::sequenceLayer(reads[[1]]$qual,
                                reads[[1]]$cigar,to = "reference")
    names(scores_aligned) <- reads[[1]]$qname

    #Fill with gaps at the end of sequences shorter than pos_end
    index_short_sequences <- which(Biostrings::width(reads_aligned) <pos_end)
    for(i in index_short_sequences){
        reads_aligned[i] <- paste0(as.character(reads_aligned[i]), paste0(rep("-", pos_end-width(reads_aligned[i])),collapse = ""),collapse = "")
        scores_aligned[i] <- paste0(as.character(scores_aligned[i]), paste0(rep("-", pos_end-width(scores_aligned[i])),collapse = ""),collapse = "")
    }
    reads_aligned <- subseq(reads_aligned, start=pos_start+1-reads[[1]]$pos,end=pos_end+1-reads[[1]]$pos)
    scores_aligned <- subseq(scores_aligned, start=pos_start+1-reads[[1]]$pos,end=pos_end+1-reads[[1]]$pos)
    # Replace deletions score's values
    if(verbose){message("Assign values to deletions\n")}
    if(verbose){pb=utils::txtProgressBar(min =0,max=length(reads_aligned),style = 3)}
    a <- str_locate_all(as.character(reads_aligned), "-+")
    for(i in seq(a)){
        if(verbose){utils::setTxtProgressBar(pb,i)}
        if(nrow(a[[i]])==0){next()}
        pos_to_replace <- IRanges(a[[i]][,"start"], a[[i]][,"end"])
        pos_to_average_start <- IRanges(a[[i]][,"start"]-1,a[[i]][,"start"]-1)
        pos_to_average_end <- IRanges(a[[i]][,"end"]+1,a[[i]][,"end"]+1)
        if(BiocGenerics::start(pos_to_average_start)[1]==0){
            BiocGenerics::start(pos_to_average_start)[1] <- BiocGenerics::end(pos_to_average_start)[1] <- BiocGenerics::start(pos_to_average_end)[1]
        }
        n <- length(pos_to_average_start)
        pos_max <- width(reads_aligned)[i]
        if(BiocGenerics::start(pos_to_average_end)[n]>pos_max){
            BiocGenerics::start(pos_to_average_end)[n] <-  BiocGenerics::start(pos_to_average_start)[n]
            BiocGenerics::end(pos_to_average_end)[n] <- BiocGenerics::start(pos_to_average_start)[n]
        }
        before <- extractAt(scores_aligned[[i]], at = pos_to_average_start)
        before <- ascii_v[as(before,"character")]
        after <- extractAt(scores_aligned[[i]], at = pos_to_average_end)
        after <- ascii_v[as(after,"character")]
        both <- cbind(before,after)
        if(gaps_weights=="mean"){
            replace_qscore <- apply(both,1, mean)
        }else if(gaps_weights=="minimum"){
            replace_qscore <- apply(both,1, min)
        }
        replace_qscore <- lapply(replace_qscore, function(x){as(x,"PhredQuality")})
        replace_qscore_v <- sapply(seq(replace_qscore), function(x){paste0(rep(replace_qscore[[x]], times=width(pos_to_replace)[x]), collapse="")})
        scores_aligned[i] <- replaceAt(scores_aligned[i],at = pos_to_replace, value = replace_qscore_v)

    }

    if(save_original_scores){
        fastq_original_scores <- QualityScaledDNAStringSet(reads_aligned,scores_aligned)
        output_original_qscores <- paste0(tools::file_path_sans_ext(output_file), "_original.fastq")
        writeQualityScaledXStringSet(fastq_original_scores, filepath = output_original_qscores)
    }
    if(verbose){message("\nCorrect QUAL values\n")}
    if(verbose){pb=utils::txtProgressBar(min =0,max=length(reads_aligned),style = 3)}
    for(i in seq(length(reads_aligned))){
        if(verbose){utils::setTxtProgressBar(pb,i)}
        mismatches <- Biostrings::compareStrings(ref_seq,reads_aligned[i])
        mismatches_positions <- stringr::str_locate_all(mismatches, pattern="\\+|\\?|\\-")
        if(nrow(mismatches_positions[[1]])==0){next()}
        mismatches_positions_Ir <- IRanges(start=mismatches_positions[[1]][,1], end=mismatches_positions[[1]][,2])
        mismatch_from <- extractAt(ref_seq,mismatches_positions_Ir)
        mismatch_to <- extractAt(reads_aligned[i],mismatches_positions_Ir)
        mismatches_q <- extractAt(scores_aligned[i],mismatches_positions_Ir)
        mismatches_df <- data.frame(nucleotide=mismatch_to[[1]],
                                    pos = mismatches_positions[[1]][,1],
                                    QUAL = unlist(as(unlist(mismatches_q),"IntegerList")),
                                    strand=unlist(reads[[1]]$strand[i]))
        mismatches_df <- mismatches_df %>%
            left_join(qtable %>% select(c("nucleotide","pos","QUAL","strand","isWT","p_SINGLe")),
                    by = c("nucleotide","pos","QUAL","strand"))
        na_single <- is.na(mismatches_df$p_SINGLe)
        if(any(na_single)){
            mismatches_df <- mismatches_df[!na_single,]
            mismatches_positions_Ir <- mismatches_positions_Ir[!na_single]
        }

        singleQUAL <- sapply(1-mismatches_df$p_SINGLe,PhredQuality)
        for(p in seq_along(mismatches_positions_Ir)){
            scores_aligned[i] <- replaceAt(scores_aligned[i],mismatches_positions_Ir[p], singleQUAL[[p]])
        }
    }
    output <- QualityScaledDNAStringSet(reads_aligned,scores_aligned)
    if(save){writeQualityScaledXStringSet(output, filepath = output_file)}
    return(output)
}


# if(F){
#     for(i in seq(length(reads_aligned))){
#         if(verbose){utils::setTxtProgressBar(pb,i)}
#         ## REPLACE DELETION SCORES
#         del_positions <- BiocGenerics::start(Biostrings::vmatchPattern("-",reads_aligned[i])[[1]])
#         ndel = length(del_positions)
#
#         if(ndel==0){next()}
#         nr = width(reads_aligned)[i]
#
#         #If they are at the beginning, I replace by the first Qscore
#         if(del_positions[1]==1){
#             n_del_start=1
#             if(length(del_positions)==1){
#                 n_del_start <- 1
#             }else{
#                 condition = del_positions[n_del_start+1]==del_positions[n_del_start]+1
#                 while(condition){
#                     n_del_start <- n_del_start+1
#                     condition1 = is.na(del_positions[n_del_start+1])
#                     if(condition1){
#                         condition=FALSE
#                     }else{
#                         condition = del_positions[n_del_start+1]==del_positions[n_del_start]+1
#                     }
#                 }
#             }
#             scores_aligned[i] <- Biostrings::replaceAt(scores_aligned[i],at=IRanges(1,n_del_start),as.character(subseq(scores_aligned[i], 1, n_del_start)))
#         }else{n_del_start=NA}
#
#         #If they are at the end, I replace by the last Qscore
#         if(del_positions[ndel]==nr){
#             #detect which are the values on the end of the string
#             Breaks <- c(0, which(diff(del_positions) != 1), length(del_positions))
#             Breaks_ini <- Breaks[length(Breaks)-1]+1
#             Breaks_end <- Breaks[length(Breaks)]
#             del_pos_ini = del_positions[Breaks_ini]
#             del_pos_end = del_positions[Breaks_end]
#             replacement <- paste0(rep(as.character(subseq(scores_aligned[i], del_pos_ini-1, del_pos_ini-1)), Breaks_end-Breaks_ini+1),collapse="")
#             scores_aligned[i] <-Biostrings::replaceAt(scores_aligned[i],
#                                                       at=IRanges(del_pos_ini,del_pos_end),
#                                                       replacement)
#             n_del_end = length(Breaks_ini:Breaks_end)
#         }else{n_del_end=NA}
#
#         if(!is.na(n_del_end)){del_positions <- del_positions[-seq(ndel-n_del_end+1,ndel)]}
#         if(!is.na(n_del_start)){del_positions <- del_positions[-seq(n_del_start)]}
#         if(length(del_positions)==0){next()}
#         ndel = length(del_positions)
#         for(j in seq_along(del_positions)){
#             jmin = j
#             while(j < ndel && del_positions[j+1]==del_positions[j]+1){ j <- j+1 }
#             jmax=j
#
#             q_upstr <- subseq(scores_aligned[i], del_positions[jmin]-1, del_positions[jmin]-1)
#             q_dwstr <- subseq(scores_aligned[i], del_positions[jmax]+1, del_positions[jmax]+1)
#
#             if(gaps_weights=="mean"){
#                 prob_upstr <- as(q_upstr,"NumericList")[[1]]
#                 prob_dwstr  <- as(q_dwstr,"NumericList")[[1]]
#                 prob <- mean(prob_dwstr,prob_upstr)
#             }else if(gaps_weights=="minimum"){
#                 prob_upstr <- as(q_upstr,"NumericList")[[1]]
#                 prob_dwstr  <- as(q_dwstr,"NumericList")[[1]]
#                 prob <- min(prob_dwstr,prob_upstr)
#             }
#
#             qscore_new <- as(prob,"PhredQuality")
#             qscore_new_vec <- paste0(rep(as.character(qscore_new), jmax-jmin+1), collapse="")
#             scores_aligned[i] <- Biostrings::replaceAt(scores_aligned[i],at=IRanges(del_positions[jmin],del_positions[jmax]),qscore_new_vec)
#         }
#
#     }
# }
