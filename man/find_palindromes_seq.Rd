\name{find_palindromes_seq}
\alias{find_palindromes_seq}
\title{Find palindrome-like inverted repeats in a DNA sequence}
\usage{
find_palindromes_seq(
  sequence,
  seq_id = NA_character_,
  min_arm = 4L,
  max_arm = 12L,
  min_loop = 0L,
  max_loop = 20L,
  max_mismatches = 0L
)
}
\arguments{
  \item{sequence}{A single DNA sequence as a character string.}
  \item{seq_id}{Optional label carried into the output.}
  \item{min_arm}{Minimum arm length to test.}
  \item{max_arm}{Maximum arm length to test.}
  \item{min_loop}{Minimum spacer/loop length to test.}
  \item{max_loop}{Maximum spacer/loop length to test.}
  \item{max_mismatches}{Maximum number of non-complementary arm positions allowed.}
}
\value{
A tibble of palindrome candidates with sequence coordinates, arm and loop sizes,
pairing statistics, and dot-bracket representation.
}
\description{
Scans a single DNA sequence for inverted repeats. Each candidate window is
translated into an RNA secondary-structure-style dot-bracket string where
matched complementary positions are shown as \code{(} and \code{)}, and mismatched
positions are shown as \code{.}.
}
\examples{
find_palindromes_seq(
  sequence = "ATGCAAGCAT",
  min_arm = 4,
  max_arm = 4,
  min_loop = 2,
  max_loop = 2,
  max_mismatches = 0
)
}
