\name{search_structure}
\alias{search_structure}
\title{Filter palindrome hits by dot-bracket pattern}
\usage{
search_structure(hits, pattern)
}
\arguments{
  \item{hits}{A tibble returned by \code{\link{find_palindromes_seq}} or \code{\link{find_palindromes_genome}}.}
  \item{pattern}{A regular expression evaluated against the \code{dot_bracket} column.}
}
\value{
A filtered tibble.
}
\description{
Uses \code{stringr::str_detect()} to retain only rows whose dot-bracket structure
matches a regular expression.
}
\examples{
hits <- find_palindromes_seq("ATGCAAGCAT", min_arm = 4, max_arm = 4, min_loop = 2, max_loop = 2)
search_structure(hits, "^\\\\({4,}\\\\.{2}\\\\){4,}$")
}
