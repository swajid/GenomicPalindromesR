\name{summarise_palindromes}
\alias{summarise_palindromes}
\title{Summarize palindrome calls}
\usage{
summarise_palindromes(hits)
}
\arguments{
  \item{hits}{A tibble returned by \code{\link{find_palindromes_seq}} or \code{\link{find_palindromes_genome}}.}
}
\value{
A tibble with counts and average structural statistics.
}
\description{
Aggregates a hit table by chromosome when present, or across all hits otherwise.
}
\examples{
hits <- find_palindromes_seq("ATGCAAGCAT", min_arm = 4, max_arm = 4, min_loop = 2, max_loop = 2)
summarise_palindromes(hits)
}
