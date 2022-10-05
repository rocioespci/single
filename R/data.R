#' Bases
#'
#' Vector A C G T -
#' @format Character vector of length 5
"bases"

#' mutation_rate
#'
#' Mutational rate matrix for error-prone PCR, obtained from GeneMorph II Random Mutagenesis Kit.
#' @format matrix size 4x5
#' @source https://www.agilent.com/cs/library/usermanuals/public/200550.pdf
"mutation_rate"

#' ASCII code
#'
#' Vector ascii
#' @format Named character vector of length 94. Names are ascii character and value is the probability of error.
"ascii_v"


#' @importFrom utils globalVariables
utils::globalVariables(c("bases", "ascii","mutation_rate","ascii_v"))
