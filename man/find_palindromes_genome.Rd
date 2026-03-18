\name{find_palindromes_genome}
\alias{find_palindromes_genome}
\title{Scan a BSgenome object for palindrome-like inverted repeats}
\usage{
find_palindromes_genome(
  genome,
  chromosomes = NULL,
  min_arm = 4L,
  max_arm = 12L,
  min_loop = 0L,
  max_loop = 20L,
  max_mismatches = 0L
)
}
\arguments{
  \item{genome}{A \code{BSgenome} object such as \code{BSgenome.Mmusculus.UCSC.mm39}.}
  \item{chromosomes}{Character vector of chromosome names. If \code{NULL}, all sequence levels are scanned.}
  \item{min_arm}{Minimum arm length to test.}
  \item{max_arm}{Maximum arm length to test.}
  \item{min_loop}{Minimum spacer/loop length to test.}
  \item{max_loop}{Maximum spacer/loop length to test.}
  \item{max_mismatches}{Maximum number of allowed mismatches per candidate.}
}
\value{
A tibble of palindrome candidates with a \code{chromosome} column.
}
\description{
Iterates across one or more chromosomes from a \code{BSgenome} object,
extracts the underlying DNA sequence, and applies \code{\link{find_palindromes_seq}}
to each chromosome.
}
\examples{
\donttest{
library(BSgenome.Mmusculus.UCSC.mm39)

find_palindromes_genome(
  genome = BSgenome.Mmusculus.UCSC.mm39,
  chromosomes = c("chr1", "chrX"),
  min_arm = 6,
  max_arm = 8,
  min_loop = 0,
  max_loop = 10,
  max_mismatches = 1
)
}
}
