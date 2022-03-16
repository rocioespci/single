#' Compute SINGLE consensus
#'
#' Main function to compute consensus after correcting reads by a SINGLE model.
#'
#' @param barcodes_file File containing the names of the reads and the barcode associated (or any grouping tag).
#' @param single_corrected_seqs Files containing the sequences corrected by SINGLE, i.e. saved by single_evaluate
#' @param column_SeqID,column_BCid Numeric. Columns where the sequences id and the barcode (or grouping tag) are, in the barcodes_file
#' @param header,dec,sep  Arguments for read.table(barcodes_file)
#' @param verbose Logical.
#' @return DNAStringSet with consensus sequences
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom utils read.table
#' @importFrom Biostrings DNAStringSet readQualityScaledDNAStringSet quality
#' @export single_consensus_byBarcode
#' @examples
#' pos_start=1
#' pos_end = 10
#' refseq_fasta = system.file("extdata", "ref_seq.fasta", package = "single")
#' ref_seq <- subseq(readDNAStringSet(refseq_fasta), pos_start,pos_end)
#' train_reads_example <- system.file("extdata", "train_seqs_500.sorted.bam",
#'                                    package = "single")
#' train <- single_train(bamfile=train_reads_example,
#'                    refseq_fasta=refseq_fasta,
#'                    rates.matrix=mutation_rate,
#'                    mean.n.mutations=5.4,
#'                    pos_start=pos_start,
#'                    pos_end=pos_end,
#'                    save_final= FALSE)
#' test_reads_example = system.file("extdata", "test_sequences.sorted.bam",
#'    package = "single")
#' corrected_reads <- single_evaluate(bamfile = test_reads_example,
#'                  single_fits = train,
#'                  ref_seq = ref_seq,
#'                  pos_start=pos_start,
#'                  pos_end=pos_end,
#'                  gaps_weights = "minimum")
#' barcodes_table_example = system.file("extdata", "Barcodes_table.txt",
#'    package = "single")
#' consensus <- single_consensus_byBarcode(
#'                  barcodes_file = barcodes_table_example,
#'                  single_corrected_seqs = corrected_reads,
#'                  verbose = FALSE)
single_consensus_byBarcode <- function(barcodes_file,
                                single_corrected_seqs,
                                column_SeqID=NULL,
                                column_BCid=NULL,
                                header=TRUE, dec=".",sep=" ", verbose=TRUE){
    bc_table <- utils::read.table(barcodes_file,header=header, dec=dec,sep=sep)
    if( (is.null(column_BCid) & !is.null(column_SeqID))){
        stop('single_consensus_byBarcode: you need to specify column_BCid')
    }
    if( (!is.null(column_BCid) & is.null(column_SeqID))){
        stop('single_consensus_byBarcode: you need to specify column_SeqID')
    }
    if(!is.null(column_BCid) | !is.null(column_SeqID)){
        bc_table <- bc_table[,c(column_SeqID,column_BCid)]
    }
    colnames(bc_table) <- c("SeqID","BCid")

    if(is.character(single_corrected_seqs)){
        single_corrected_seqs <- readQualityScaledDNAStringSet(single_corrected_seqs)
    }
    barcodes = unique(bc_table$BCid)

    #Compute consensus for each barcode using all sequences available and using the three available methods
    consensus_sequences <- DNAStringSet()
    if(verbose){ pb <- utils::txtProgressBar(0,length(barcodes),style=3) }
    for(bc in seq_along(barcodes)){
        bc_table_aux <- bc_table %>%
            filter(BCid==barcodes[bc])
        index <- which(names(single_corrected_seqs) %in% bc_table_aux$SeqID)
        seqs_in_barcode <- single_corrected_seqs[index]
        aux_seqs <- sapply(as.character(seqs_in_barcode),strsplit, split="")
        aux_quals <- 1-as(quality(seqs_in_barcode),"NumericList")
        aux_pos <- lapply(aux_seqs, seq_along)
        data_barcode = data.frame(nucleotide = unlist(aux_seqs),
                                    p_SINGLe=unlist(aux_quals),
                                    pos=unlist(aux_pos) )
        rownames(data_barcode) <- NULL
        consensus_sequences[[bc]] <- weighted_consensus(df = data_barcode, cutoff_prob = 0)
        if(verbose){utils::setTxtProgressBar(pb,bc)}
    }
    return(consensus_sequences)
}

