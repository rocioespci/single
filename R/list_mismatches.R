#' Lists mismatches between two DNAstrings
#'
#' This is an auxiliary function in single package, to list the mutations of two DNAstrings.
#' @param ref DNAString, reference sequence.
#' @param seq DNAString, target sequence, same length as ref.
#' @return Character vector containing Nucleotide in ref Position Nucleotide in seq. If ref and seq are equal, it returns NA.
#' @importFrom stringr str_locate_all str_sub
#' @importFrom Biostrings compareStrings DNAString
#' @export list_mismatches
#' @examples
#' ref = Biostrings::DNAString("AAAA")
#' seq = Biostrings::DNAString("AGAT")
#' list_mismatches(ref,seq)
#' list_mismatches(ref,ref)
list_mismatches <- function(ref,seq){
    compare <- Biostrings::compareStrings(ref,seq)
    pos     <- stringr::str_locate_all(compare,"[\\?\\-\\+]")
    if(length(pos[[1]])==0){y=NA}else{
        from    <- sapply(pos,function(x){stringr::str_sub(as.character(ref),x)})
        to      <- sapply(pos,function(x){stringr::str_sub(as.character(seq),x)})
        y       <- paste0(from,pos[[1]][,1],to)
    }
    return(y)
}
