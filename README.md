# GenomicPalindromesR

`GenomicPalindromesR` finds DNA palindromes and inverted repeats in mouse genomic
sequence by converting candidate windows into an RNA secondary-structure-style
dot-bracket representation such as `((((....))))` and then searching those
strings with tidyverse tools like `dplyr::filter()` and `stringr::str_detect()`.

## Installation

```r
# install.packages("pak")
pak::pak(c(
  "tidyverse",
  "BiocManager"
))

BiocManager::install(c(
  "BSgenome",
  "Biostrings",
  "GenomeInfoDb",
  "BSgenome.Mmusculus.UCSC.mm39"
))

# from a local clone
pak::local_install(".")
```

## Why dot-bracket for DNA palindromes?

A perfect inverted repeat can be represented the same way an RNA hairpin is:

- left arm: `((((`
- loop/spacer: `....`
- right arm: `))))`

A near-palindrome with broken complementarity becomes something like
`((.. ... ..))`, where mismatched positions are marked as dots instead of
parentheses. This lets you use regular expressions on structure alone.

## Example: scan a sequence

```r
library(GenomicPalindromesR)
library(dplyr)

hits <- find_palindromes_seq(
  sequence = "ATGCAAGCAT",
  min_arm = 4,
  max_arm = 4,
  min_loop = 2,
  max_loop = 2,
  max_mismatches = 0
)

hits %>%
  search_structure("^\\({4}\\.{2}\\){4}$")
```

## Example: scan the mouse genome

```r
library(GenomicPalindromesR)
library(BSgenome.Mmusculus.UCSC.mm39)
library(dplyr)

mm39_hits <- find_palindromes_genome(
  genome = BSgenome.Mmusculus.UCSC.mm39,
  chromosomes = c("chr1", "chr7", "chrX"),
  min_arm = 6,
  max_arm = 12,
  min_loop = 0,
  max_loop = 20,
  max_mismatches = 1
)

mm39_hits %>%
  search_structure("^\\({6,}\\.{0,20}\\){6,}$") %>%
  summarise_palindromes()
```

## Main functions

- `find_palindromes_seq()` scans a single DNA sequence.
- `find_palindromes_genome()` scans chromosomes from a `BSgenome` object.
- `search_structure()` filters hits using a regex over the dot-bracket string.
- `summarise_palindromes()` aggregates hits by chromosome and structure class.

## Returned columns

Typical output includes:

- `chromosome`
- `start`
- `end`
- `width`
- `arm_length`
- `loop_length`
- `matched_pairs`
- `mismatches`
- `pairing_fraction`
- `dot_bracket`
- `left_arm`
- `loop_sequence`
- `right_arm`
- `window_sequence`

## Notes

This package intentionally emphasizes readability and tidyverse workflows over
maximal genome-scale optimization. For whole-genome scans, start with a subset
of chromosomes or narrower arm/loop ranges.


## Additional Notes
The original algorithm of this, essentially, the idea that a whole genome should be converted into RNA structure dot parenthesis notation (through software that does that, written in C++), (not RNA but DNA!) was and is an original idea by me because I worked on RNA structure dot parenthesis notation during a sub-project in my iGEM team (2011) with my Bioinformatics Prof (now a Group Leader at BMS). My case being, "why write a palindrome solver from scratch when one has already been written." And it worked. (Also fun fact: I originally wrote a palindrome solver from scratch in junior year of high school in Java. Hence my approach: been there, done that!)

# Pull Request Inquiries

Open to continue to develop this R Package / Publication in the near future. Please do kindly make a pull request!

Teşekkürler <-> Thank You

Although, on second thought, I've moved SO BEYOND THIS that it really doesn't matter, umm yeah ... still onward to more interesting things.

We don't care about this, it's done with. ARCHIVED.

![img](rewrite-it-in-rust-meme.png)

# Press Inquiries
Sure, please email me @ s a n a w g s @ g m a i l . c o m
