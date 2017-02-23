context("nltt_diff_exact_norm_brts")


test_that("nltt_diff_exact_norm_brts: use", {

  if (1 == 2) {
    set.seed(42)
    p <- ape::rcoal(5)
    p$root.edge <- 0.1
    q <- ape::rcoal(5)
    p$root.edge <- 0.2
    nLTT::nLTTstat_exact(p, q, ignore_stem = FALSE)
    nLTT::nLTTstat_exact(p, q, ignore_stem = TRUE)
  }
  # High warning level
  options(warn = 2)

  if (2 == 3) {
    n <- 10
    b_times_N <- (seq(1, n) / n)
    lineages_N <- b_times_N
    b_times2_N <- b_times_N * b_times_N
    lineages2_N <- b_times2_N

    measured <- nLTT::nltt_diff_exact_norm_brts(
      b_times_N = b_times_N,
      lineages_N = lineages_N,
      b_times2_N = b_times2_N,
      lineages2_N = lineages2_N,
      distance_method = "abs"
    )
  }
})
