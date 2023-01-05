#' Compute SINGLE consensus
#'
#' Main function to compute consensus after correcting reads by a SINGLE model.
#'
#' @param barcodes_table data.frame or file name containing the names of the reads and the barcode associated (or any grouping tag).
#' @param sequences QualityScaledDNAStringSet or fastq file name. Contains sequences from which compute weighted consensus.
#' @param readID_col,bcID_col Numeric. Columns where the reads id and the barcode (or grouping tag) are, in the barcodes_table
#' @param header,dec,sep  Arguments for read.table(barcodes_table)
#' @param verbose Logical.
#' @return DNAStringSet with consensus sequences
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom utils read.table
#' @importFrom methods as
#' @importFrom Biostrings DNAStringSet readQualityScaledDNAStringSet quality readDNAStringSet
#' @export single_consensus_byBarcode
#' @examples
#' pos_start=1
#' pos_end = 100
#' barcodes_file =  system.file("extdata", "Barcodes_table.txt",package = "single")
#' reads_single = system.file("extdata", "corrected_seqs.fastq", package = "single")
#' single_consensus_byBarcode(barcodes_file,reads_single, verbose = FALSE)
single_consensus_byBarcode <- function(barcodes_table,sequences,
                                readID_col=1,bcID_col=2,
                                header=TRUE, dec=".",sep=" ", verbose=TRUE){
    if(is.character(barcodes_table)){
        barcodes_table <- utils::read.table(barcodes_table,
                                     header=header,dec=dec,sep=sep)
    }
    barcodes_table  <- barcodes_table[,c(readID_col,bcID_col)]
    colnames(barcodes_table) <- c("readID","bcID")

    if(is.character(sequences)){
        sequences <- readQualityScaledDNAStringSet(sequences)
    }
    #Intersect sequences both in sequences and barcodes_table
    intersection_names <- intersect(barcodes_table$readID,names(sequences))
    barcodes_table <- barcodes_table %>% filter(readID %in% intersection_names)
    sequences <- sequences[intersection_names]

    #Compute consensus for each barcode
    barcodes <- unique(barcodes_table$bcID)
    consensus_sequences <- DNAStringSet()
    names_barcodes <- rep(NA,length(barcodes))
    if(verbose){ pb <- utils::txtProgressBar(0,length(barcodes),style=3) }
    for(bc in seq_along(barcodes)){
        bc_table_aux <- barcodes_table %>%
            filter(bcID==barcodes[bc])
        if(nrow(bc_table_aux)==1){
            consensus_sequences[[bc]] <- sequences[bc_table_aux$readID][[1]]
        }
        seqs_in_barcode <- sequences[bc_table_aux$readID]
        aux_seqs <- sapply(as.character(seqs_in_barcode),strsplit, split="")
        aux_quals <- 1-as(quality(seqs_in_barcode),"NumericList")
        aux_pos <- lapply(aux_seqs, seq_along)
        data_barcode = data.frame(nucleotide = unlist(aux_seqs),
                                    p_SINGLe=unlist(aux_quals),
                                    pos=unlist(aux_pos) )
        rownames(data_barcode) <- NULL
        consensus_sequences[[bc]] <- weighted_consensus(df = data_barcode, cutoff_prob = 0)
        names_barcodes[bc] <- barcodes[bc]
        if(verbose){utils::setTxtProgressBar(pb,bc)}
    }
    names(consensus_sequences) <- names_barcodes
    return(consensus_sequences)
}

