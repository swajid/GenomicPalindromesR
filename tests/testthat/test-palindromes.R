test_that("perfect palindrome receives fully bracketed structure", {
  hits <- find_palindromes_seq(
    sequence = "ATGCAT",
    min_arm = 3,
    max_arm = 3,
    min_loop = 0,
    max_loop = 0,
    max_mismatches = 0
  )

  expect_equal(nrow(hits), 1L)
  expect_equal(hits$dot_bracket[[1]], "((()))")
  expect_equal(hits$mismatches[[1]], 0)
  expect_equal(hits$matched_pairs[[1]], 3)
})

test_that("mismatched positions become dots in dot-bracket output", {
  hits <- find_palindromes_seq(
    sequence = "ATGCAA",
    min_arm = 3,
    max_arm = 3,
    min_loop = 0,
    max_loop = 0,
    max_mismatches = 3
  )

  expect_equal(nrow(hits), 1L)
  expect_true(grepl("\\.", hits$dot_bracket[[1]]))
  expect_equal(hits$mismatches[[1]], 1)
})

test_that("structure regex filters work", {
  hits <- find_palindromes_seq(
    sequence = "ATGCAAGCAT",
    min_arm = 4,
    max_arm = 4,
    min_loop = 2,
    max_loop = 2,
    max_mismatches = 0
  )

  filtered <- search_structure(hits, "^\\({4}\\.{2}\\){4}$")
  expect_equal(nrow(filtered), 1L)
})

test_that("summary preserves chromosome grouping when present", {
  x <- tibble::tibble(
    chromosome = c("chr1", "chr1", "chr2"),
    arm_length = c(4L, 4L, 6L),
    loop_length = c(2L, 2L, 0L),
    matched_pairs = c(4L, 3L, 6L),
    mismatches = c(0L, 1L, 0L),
    pairing_fraction = c(1, 0.75, 1)
  )

  out <- summarise_palindromes(x)
  expect_true(all(c("chromosome", "n_hits") %in% names(out)))
})
