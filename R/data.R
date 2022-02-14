#' Bases
#'
#' Vector A C G T -
#' @format Character vector of length 5
"bases"

#' ascii
#'
#' Equivalences of ascii code. Columns are Number, (ascii) Symbol, Q (score) and P (probability of error)
#' @format data.frame size 94x4
"ascii"

#' mutation_rate
#'
#' Mutational rate matrix for error-prone PCR, obtained from GeneMorph II Random Mutagenesis Kit.
#' @format matrix size 4x5
#' @source https://www.agilent.com/cs/library/usermanuals/public/200550.pdf
"mutation_rate"

#' @importFrom utils globalVariables
utils::globalVariables(c("bases", "ascii","mutation_rate"))

