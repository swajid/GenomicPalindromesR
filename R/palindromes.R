#' Internal helpers for palindrome scanning
#'
#' These helpers normalize DNA input, compute reverse-complement-aware pairing,
#' and construct RNA-style dot-bracket strings.
#' @keywords internal
NULL

.normalize_sequence <- function(sequence) {
  if (length(sequence) != 1L || is.na(sequence)) {
    rlang::abort("`sequence` must be a single non-missing character string.")
  }

  sequence <- toupper(gsub("\\s+", "", sequence))

  if (!nzchar(sequence)) {
    rlang::abort("`sequence` must not be empty.")
  }

  if (grepl("[^ACGTU]", sequence)) {
    rlang::abort("`sequence` may contain only A, C, G, T, or U.")
  }

  chartr("U", "T", sequence)
}

.base_complement <- function(x) {
  lookup <- c(A = "T", T = "A", C = "G", G = "C")
  unname(lookup[x])
}

.split_bases <- function(sequence) {
  strsplit(sequence, split = "", fixed = TRUE)[[1]]
}

.build_dot_bracket <- function(left_arm, loop_sequence, right_arm) {
  left_bases <- .split_bases(left_arm)
  right_bases <- .split_bases(right_arm)

  if (length(left_bases) != length(right_bases)) {
    rlang::abort("`left_arm` and `right_arm` must have the same length.")
  }

  right_partners <- rev(right_bases)
  matches <- .base_complement(left_bases) == right_partners

  left_marks <- ifelse(matches, "(", ".")
  right_marks <- rev(ifelse(matches, ")", "."))

  list(
    matched_pairs = sum(matches),
    mismatches = sum(!matches),
    pairing_fraction = mean(matches),
    dot_bracket = paste0(
      paste0(left_marks, collapse = ""),
      stringr::str_dup(".", nchar(loop_sequence)),
      paste0(right_marks, collapse = "")
    )
  )
}

#' Find palindrome-like inverted repeats in a DNA sequence
#'
#' Scans a single DNA sequence for inverted repeats. Each candidate window is
#' translated into an RNA secondary-structure-style dot-bracket string where
#' matched complementary positions are shown as `(` and `)`, and mismatched
#' positions are shown as `.`.
#'
#' @param sequence A single DNA sequence as a character string.
#' @param seq_id Optional label carried into the output.
#' @param min_arm Minimum arm length to test.
#' @param max_arm Maximum arm length to test.
#' @param min_loop Minimum spacer/loop length to test.
#' @param max_loop Maximum spacer/loop length to test.
#' @param max_mismatches Maximum number of non-complementary arm positions
#'   allowed.
#'
#' @return A tibble of palindrome candidates with sequence coordinates, arm and
#'   loop sizes, pairing statistics, and dot-bracket representation.
#' @export
#'
#' @examples
#' find_palindromes_seq(
#'   sequence = "ATGCAAGCAT",
#'   min_arm = 4,
#'   max_arm = 4,
#'   min_loop = 2,
#'   max_loop = 2,
#'   max_mismatches = 0
#' )
find_palindromes_seq <- function(sequence,
                                 seq_id = NA_character_,
                                 min_arm = 4L,
                                 max_arm = 12L,
                                 min_loop = 0L,
                                 max_loop = 20L,
                                 max_mismatches = 0L) {
  sequence <- .normalize_sequence(sequence)

  min_arm <- as.integer(min_arm)
  max_arm <- as.integer(max_arm)
  min_loop <- as.integer(min_loop)
  max_loop <- as.integer(max_loop)
  max_mismatches <- as.integer(max_mismatches)

  if (min_arm < 1L || max_arm < min_arm) {
    rlang::abort("Arm lengths must satisfy `1 <= min_arm <= max_arm`.")
  }

  if (min_loop < 0L || max_loop < min_loop) {
    rlang::abort("Loop lengths must satisfy `0 <= min_loop <= max_loop`.")
  }

  if (max_mismatches < 0L) {
    rlang::abort("`max_mismatches` must be non-negative.")
  }

  n_bases <- nchar(sequence)
  results <- list()
  idx <- 1L

  for (arm_length in seq.int(min_arm, max_arm)) {
    for (loop_length in seq.int(min_loop, max_loop)) {
      width <- 2L * arm_length + loop_length

      if (width > n_bases) {
        next
      }

      starts <- seq_len(n_bases - width + 1L)

      for (start in starts) {
        end <- start + width - 1L
        left_start <- start
        left_end <- start + arm_length - 1L
        loop_start <- left_end + 1L
        loop_end <- loop_start + loop_length - 1L
        right_start <- loop_end + 1L
        right_end <- end

        left_arm <- substr(sequence, left_start, left_end)
        loop_sequence <- if (loop_length == 0L) "" else substr(sequence, loop_start, loop_end)
        right_arm <- substr(sequence, right_start, right_end)

        structure <- .build_dot_bracket(
          left_arm = left_arm,
          loop_sequence = loop_sequence,
          right_arm = right_arm
        )

        if (structure$mismatches <= max_mismatches) {
          results[[idx]] <- tibble::tibble(
            seq_id = seq_id,
            start = start,
            end = end,
            width = width,
            arm_length = arm_length,
            loop_length = loop_length,
            matched_pairs = structure$matched_pairs,
            mismatches = structure$mismatches,
            pairing_fraction = structure$pairing_fraction,
            dot_bracket = structure$dot_bracket,
            left_arm = left_arm,
            loop_sequence = loop_sequence,
            right_arm = right_arm,
            window_sequence = substr(sequence, start, end)
          )
          idx <- idx + 1L
        }
      }
    }
  }

  if (length(results) == 0L) {
    return(
      tibble::tibble(
        seq_id = character(),
        start = integer(),
        end = integer(),
        width = integer(),
        arm_length = integer(),
        loop_length = integer(),
        matched_pairs = integer(),
        mismatches = integer(),
        pairing_fraction = numeric(),
        dot_bracket = character(),
        left_arm = character(),
        loop_sequence = character(),
        right_arm = character(),
        window_sequence = character()
      )
    )
  }

  dplyr::bind_rows(results) |>
    dplyr::arrange(
      dplyr::desc(.data$pairing_fraction),
      dplyr::desc(.data$arm_length),
      .data$loop_length,
      .data$start
    )
}

#' Scan a BSgenome object for palindrome-like inverted repeats
#'
#' Iterates across one or more chromosomes from a `BSgenome` object, extracts
#' the underlying DNA sequence, and applies [find_palindromes_seq()] to each
#' chromosome.
#'
#' @param genome A `BSgenome` object such as
#'   `BSgenome.Mmusculus.UCSC.mm39`.
#' @param chromosomes Character vector of chromosome names. If `NULL`, all
#'   sequence levels in `genome` are scanned.
#' @param min_arm Minimum arm length to test.
#' @param max_arm Maximum arm length to test.
#' @param min_loop Minimum spacer/loop length to test.
#' @param max_loop Maximum spacer/loop length to test.
#' @param max_mismatches Maximum number of allowed mismatches per candidate.
#'
#' @return A tibble of palindrome candidates with a `chromosome` column.
#' @export
#'
#' @examples
#' \dontrun{
#' library(BSgenome.Mmusculus.UCSC.mm39)
#'
#' find_palindromes_genome(
#'   genome = BSgenome.Mmusculus.UCSC.mm39,
#'   chromosomes = c("chr1", "chrX"),
#'   min_arm = 6,
#'   max_arm = 8,
#'   min_loop = 0,
#'   max_loop = 10,
#'   max_mismatches = 1
#' )
#' }
find_palindromes_genome <- function(genome,
                                    chromosomes = NULL,
                                    min_arm = 4L,
                                    max_arm = 12L,
                                    min_loop = 0L,
                                    max_loop = 20L,
                                    max_mismatches = 0L) {
  if (!inherits(genome, "BSgenome")) {
    rlang::abort("`genome` must inherit from class 'BSgenome'.")
  }

  if (is.null(chromosomes)) {
    chromosomes <- GenomeInfoDb::seqlevels(genome)
  }

  chromosomes <- as.character(chromosomes)

  missing_chroms <- setdiff(chromosomes, GenomeInfoDb::seqlevels(genome))
  if (length(missing_chroms) > 0L) {
    rlang::abort(
      paste0(
        "These chromosomes are not present in `genome`: ",
        paste(missing_chroms, collapse = ", ")
      )
    )
  }

  purrr::map_dfr(chromosomes, function(chromosome) {
    seq_chr <- as.character(BSgenome::getSeq(genome, names = chromosome))

    find_palindromes_seq(
      sequence = seq_chr,
      seq_id = chromosome,
      min_arm = min_arm,
      max_arm = max_arm,
      min_loop = min_loop,
      max_loop = max_loop,
      max_mismatches = max_mismatches
    ) |>
      dplyr::mutate(chromosome = chromosome, .before = 1)
  })
}

#' Filter palindrome hits by dot-bracket pattern
#'
#' Uses `stringr::str_detect()` to retain only rows whose dot-bracket structure
#' matches a regular expression.
#'
#' @param hits A tibble returned by [find_palindromes_seq()] or
#'   [find_palindromes_genome()].
#' @param pattern A regular expression evaluated against the `dot_bracket`
#'   column.
#'
#' @return A filtered tibble.
#' @export
#'
#' @examples
#' hits <- find_palindromes_seq("ATGCAAGCAT", min_arm = 4, max_arm = 4, min_loop = 2, max_loop = 2)
#' search_structure(hits, "^\\({4,}\\.{2}\\){4,}$")
search_structure <- function(hits, pattern) {
  if (!"dot_bracket" %in% names(hits)) {
    rlang::abort("`hits` must contain a `dot_bracket` column.")
  }

  dplyr::filter(hits, stringr::str_detect(.data$dot_bracket, pattern))
}

#' Summarize palindrome calls
#'
#' Aggregates a hit table by chromosome when present, or across all hits
#' otherwise.
#'
#' @param hits A tibble returned by [find_palindromes_seq()] or
#'   [find_palindromes_genome()].
#'
#' @return A tibble with counts and average structural statistics.
#' @export
#'
#' @examples
#' hits <- find_palindromes_seq("ATGCAAGCAT", min_arm = 4, max_arm = 4, min_loop = 2, max_loop = 2)
#' summarise_palindromes(hits)
summarise_palindromes <- function(hits) {
  required_cols <- c("arm_length", "loop_length", "matched_pairs", "mismatches", "pairing_fraction")
  missing_cols <- setdiff(required_cols, names(hits))

  if (length(missing_cols) > 0L) {
    rlang::abort(
      paste0(
        "`hits` is missing required columns: ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }

  if ("chromosome" %in% names(hits)) {
    hits |>
      dplyr::group_by(.data$chromosome, .data$arm_length, .data$loop_length) |>
      dplyr::summarise(
        n_hits = dplyr::n(),
        mean_matched_pairs = mean(.data$matched_pairs),
        mean_mismatches = mean(.data$mismatches),
        mean_pairing_fraction = mean(.data$pairing_fraction),
        .groups = "drop"
      ) |>
      dplyr::arrange(
        dplyr::desc(.data$n_hits),
        dplyr::desc(.data$mean_pairing_fraction),
        dplyr::desc(.data$arm_length)
      )
  } else {
    hits |>
      dplyr::group_by(.data$arm_length, .data$loop_length) |>
      dplyr::summarise(
        n_hits = dplyr::n(),
        mean_matched_pairs = mean(.data$matched_pairs),
        mean_mismatches = mean(.data$mismatches),
        mean_pairing_fraction = mean(.data$pairing_fraction),
        .groups = "drop"
      ) |>
      dplyr::arrange(
        dplyr::desc(.data$n_hits),
        dplyr::desc(.data$mean_pairing_fraction),
        dplyr::desc(.data$arm_length)
      )
  }
}
