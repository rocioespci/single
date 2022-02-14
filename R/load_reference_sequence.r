#' Load reference sequence from fasta file
#'
#' @param file fasta file with reference sequence
#' @return The reference sequence as a vector of characters
#' @examples
#' load_reference_sequence(system.file("extdata", "ref_seq.fasta", package = "single"))
#' @export load_reference_sequence
load_reference_sequence <- function(file){
    reference_sequence = readLines(file)

    names_lines = which(substr(reference_sequence,1,1)==">")
    if(length(names_lines)>1){stop('Reference file has more than one sequence')}
    if(length(names_lines)==0){stop('Reference file does not contain a sequence in fasta format')}

    reference_sequence = reference_sequence[-names_lines]
    if(length(reference_sequence)>1){
        reference_sequence = paste0(reference_sequence, collapse = "")
    }

    reference_sequence <- toupper(reference_sequence)
    reference_sequence <- strsplit(reference_sequence, split="")[[1]]

    return(reference_sequence)
}
