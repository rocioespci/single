\name{parse_labels}
\alias{parse_labels}
\title{parse_labels}
\usage{
parse_labels(strings, split, element,name=F)
}
\description{
Parses a vector of strings (strsplit plus keep one element.)
}
\arguments{
  \item{strings}{Strings vector, to be parsed.}
  \item{split}{Character. Split strings by this character.}
  \item{element}{Numeric. Which element to keep after split.}
  \item{name}{Logical. Should it return a named vector?}
}
\value{
Returns a character vector. At each element of strings, it splits it by 'split' and keeps element 'element' of the splitted string.
}
\examples{
x = c("h_1","h_2","h_3")
parse_labels(x,split="_", element=2, name=F)
}
