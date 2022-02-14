#' Load reference sequence from fasta file
#'
#' @param file fasta file with reference sequence
#' @return The reference sequence as a vector of characters
#' @examples
#' load_ref_seq(system.file("extdata", "ref_seq.fasta", package = "single"))
#' @export load_ref_seq
load_ref_seq <- function(file){
    ref_seq = readLines(file)

    names_lines = which(substr(ref_seq,1,1)==">")
    if(length(names_lines)>1){
        stop('Reference file has more than one sequence')
    }
    if(length(names_lines)==0){
        stop('Reference file does not contain a sequence in fasta format')
    }

    ref_seq = ref_seq[-names_lines]
    if(length(ref_seq)>1){
        ref_seq = paste0(ref_seq, collapse = "")
    }

    ref_seq <- toupper(ref_seq)
    ref_seq <- strsplit(ref_seq, split="")[[1]]

    return(ref_seq)
}
