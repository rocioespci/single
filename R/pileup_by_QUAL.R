#' Pileup by QUAL
#'
#' To explain
#' @param bam_file Bam file to pile up
#' @param QUAL_values Numeric vector. QUAL values to analyze in the data.
#' @param pos_start Numeric. Position to start analyzing, counting starts from 1 and it refers to reference used for minimap2 alignment.
#' @param pos_end  Numeric. Position to stop analyzing, counting starts from 1 and it refers to reference used for minimap2 alignment.
#' @return data.frame with columns strand,pos,nucleotide,QUAL,countss
#' @importFrom Rsamtools  PileupParam ScanBamParam pileup
#' @importFrom GenomicAlignments seqnames strand
#' @importFrom IRanges pos
#' @export pileup_by_QUAL
#' @examples
#' pos_start <- 1
#' pos_end <- 10
#' refseq_fasta <- system.file("extdata", "ref_seq.fasta", package = "single")
#' train_reads_example <- system.file("extdata", "train_seqs_500.sorted.bam",
#'                                    package = "single")
#' counts_pnq <- pileup_by_QUAL(bam_file=train_reads_example,
#'     pos_start=pos_start,pos_end=pos_end)
#' head(counts_pnq)
pileup_by_QUAL <- function(bam_file, QUAL_values=seq(93,0),pos_start=NA,pos_end=NA){
    QUAL_values <- sort(QUAL_values, decreasing = TRUE)
    cond = FALSE
    for(q in QUAL_values){
        p_param <- Rsamtools::PileupParam(max_depth=1000000,min_base_quality=q,
                               distinguish_strands=TRUE,
                               distinguish_nucleotides=TRUE,
                               include_deletions = TRUE,
                               min_mapq = 0,
                               min_nucleotide_depth=0,
                               min_minor_allele_depth=0)
        s_param <- Rsamtools::ScanBamParam(mapqFilter=0)
        pileup_larger_than_q <- Rsamtools::pileup(BamFile(bam_file),
                                       scanBamParam = s_param,
                                       pileupParam = p_param )%>%
                                select(-seqnames)
        if(!is.na(pos_end)){
            pileup_larger_than_q <- pileup_larger_than_q %>%
                filter(pos <=pos_end )
        }
        if(!is.na(pos_start)){
            pileup_larger_than_q <- pileup_larger_than_q %>%
                filter(pos>=pos_start) %>%
                mutate(pos = pos-pos_start+1)
        }
        if(nrow(pileup_larger_than_q)==0){next()}
        if(cond){
            pileup_q <- left_join(x=pileup_larger_than_q,y=pileup_memory,
                                  by= c("pos", "strand", "nucleotide")) %>%
                mutate(count.y=ifelse(is.na(count.y),0,count.y)) %>%
                mutate(count=  count.x - count.y) %>%
                arrange(pos) %>%
                select(pos,strand,nucleotide,count)%>%
                mutate(QUAL=q) %>%
                filter(count!=0)

            pileup_QUAL <- rbind(pileup_QUAL,pileup_q)
            pileup_memory <- pileup_larger_than_q

        }else{
            pileup_memory <- pileup_larger_than_q
            pileup_QUAL <- pileup_larger_than_q  %>%
                mutate(QUAL=q)
        }
        cond=TRUE
    }
    pileup_QUAL <- pileup_QUAL %>%
        arrange(strand,pos,nucleotide,QUAL) %>%
        mutate(strand=droplevels(strand),nucleotide=droplevels(nucleotide))
    return(pileup_QUAL)
}
