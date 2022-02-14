#' Compute SINGLE consensus
#'
#' Main function to compute consensus after correcting reads by a SINGLE model.
#'
#' @param barcodes_table_file File containing the names of the reads and the barcode associated (or any grouping tag).
#' @param reference_sequence_fastafile Fasta file containing reference sequence
#' @param single_corrected_files Files containing the sequences corrected by SINGLE, i.e. saved by single_evaluate_model.
#' @param file_out_consensus String. File to store results.
#' @param column_id,column_barcode Numeric. Columns where the sequences id and the barcode (or grouping tag) are, in the barcodes_table_file
#' @param header,dec,sep  Arguments for read.table(barcodes_table_file)
#' @param verbose Logical.
#' @return This function returns 0. Results are saved in file_out_consensus
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom utils read.table
#' @export single_consensus_byBarcode
#' @examples
#' barcodes_table_file = system.file("extdata", "example_barcodes_table.txt", package = "single")
#' ref_seq_file = system.file("extdata", "ref_seq.fasta", package = "single")
#' single_corrected_files = system.file("extdata", "example_corrected.txt", package="single")
#' file_out = tempfile("example_consensus.txt")
#' single_consensus_byBarcode(barcodes_table_file,ref_seq_file,single_corrected_files,file_out)
single_consensus_byBarcode <- function(barcodes_table_file,reference_sequence_fastafile,single_corrected_files,file_out_consensus,
                                       column_id=NULL,column_barcode=NULL,
                                       header=TRUE, dec=".",sep=" ", verbose=TRUE){
    barcodes_table <- utils::read.table(barcodes_table_file,header=header, dec=dec,sep=sep)
    if( (is.null(column_barcode) & !is.null(column_id))){stop('single_consensus_byBarcode: you need to specify column_barcode')}
    if( (!is.null(column_barcode) & is.null(column_id))){stop('single_consensus_byBarcode: you need to specify column_id')}

    if(!is.null(column_barcode) | !is.null(column_id)){
        barcodes_table <- barcodes_table[,c(column_id,column_barcode)]
    }
    colnames(barcodes_table) <- c("SeqID","BCsequence")

    reference <- load_reference_sequence(reference_sequence_fastafile)
    reads_corrected <- lapply(single_corrected_files,read.table, header=TRUE)
    reads_corrected <- do.call(rbind,reads_corrected)
    reads_corrected <- dplyr::left_join(reads_corrected, barcodes_table %>% select(.data$SeqID,.data$BCsequence) , by="SeqID") %>%
    dplyr::filter(!is.na(.data$BCsequence))
    reads_corrected$p_right_priors_model[is.na(reads_corrected$p_right_priors_model)]<-1

    rm(barcodes_table)
    barcodes = unique(reads_corrected$BCsequence)

    #Compute consensus for each barcode using all sequences available and using the three available methods
    con_out <- file(file_out_consensus,"w")
    cat("Barcode\tmutations_by_single\n", file=con_out)
    if(verbose){ pb <- utils::txtProgressBar(0,length(barcodes),style=3) }
    for(bc in seq_along(barcodes)){
        data_bc <- reads_corrected %>%
            dplyr::filter(.data$BCsequence==barcodes[bc])
        cons_single <- compute_consensus_sequence(df = data_bc[,c("nucleotide","p_right_priors_model","position")], cutoff_prob = 0)
        muts_single <- c(detect_mutations(sequence =  cons_single, reference = reference))
        cat(as.character(barcodes[bc]),":",muts_single, "\n",file=con_out)
        if(verbose){utils::setTxtProgressBar(pb,bc)}
    }
    close(con_out)
    return(0)
}

