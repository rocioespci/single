#' Detect mutations beween sequence and reference
#'
#' This is an auxiliary function in single package. It Compares the vectors sequence and reference, and returns the differences.
#' @param sequence Vector of characters, sequence you want to compare to reference
#' @param reference Vector of characters, reference sequence.
#' @return Character vector, in each site the difference between both vectors in the form reference - index - sequence
#' @export detect_mutations
#' @examples
#' detect_mutations(c("A","A","A","A"),c("B","A","B","A"))
detect_mutations <- function(sequence,reference){
      ind <- which(sequence != reference)
      muts <- paste0(reference[ind],ind,sequence[ind])
      return(muts)
}
