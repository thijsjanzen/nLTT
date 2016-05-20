context("stretch_nltt_matrix")

test_that("stretch_nltt_matrix #1", {
  # t   N      t   N
  # 0.0 0.5    0.0 0.5
  #         -> 0.5 0.5
  # 1.0 1.0    1.0 1.0
  m <- matrix(c(c(0.0, 1.0), c(0.5, 1.0)), ncol = 2, nrow = 2)
  colnames(m) <- c("t", "N")
  expected <- matrix(
    c(
      seq(0.0, 1.0, 0.5),
      c(0.5, 0.5, 1.0)
    ),
    ncol = 2, nrow = 3
  )
  result <- stretch_nltt_matrix(m = m, dt = 0.5, step_type = "lower")
  if (!identical(result, expected)) {
    print("ERROR")
    print("result:")
    print(result)
    print("expected:")
    print(expected)
  }
  expect_equal(identical(result, expected), TRUE)
})

test_that("stretch_nltt_matrix #2", {
  # t   N      t   N
  # 0.0 0.3    0.00 0.3
  # 0.4 0.5    0.25 0.3
  #            0.50 0.5
  #         -> 0.75 0.5
  # 1.0 1.0    1.00 1.0
  m <- matrix(c(c(0.0, 0.4, 1.0), c(0.3, 0.5, 1.0)), ncol = 2, nrow = 3)
  colnames(m) <- c("t", "N")
  expected <- matrix(
    c(
      seq(0.0, 1.0, 0.25),
      c(0.3, 0.3, 0.5, 0.5, 1.0)
    ),
    ncol = 2, nrow = 5
  )
  result <- stretch_nltt_matrix(m = m, dt = 0.25, step_type = "lower")
  if (!identical(result, expected)) {
    print("ERROR")
    print("result:")
    print(result)
    print("expected:")
    print(expected)
  }
  expect_equal(identical(result, expected), TRUE)
})






test_that("stretch_nltt_matrix from vignette, ", {
  # t   N      t   N
  # 0.0 0.3    0.00 0.3
  # 0.4 0.5    0.25 0.3
  #            0.50 0.5
  #         -> 0.75 0.5
  # 1.0 1.0    1.00 1.0
  newick <- "((A:1,B:1):1,(C:1,D:1):1);"
  phylogeny <- ape::read.tree(text = newick)
  nltt <- nLTT::get_phylogeny_nltt_matrix(phylogeny)
  result <- stretch_nltt_matrix(nltt, dt = 0.25, step_type = "upper")

  expected <- matrix(
    c(
      c(0.00, 0.25, 0.50, 0.75, 1.00),
      c(0.50, 0.50, 1.00, 1.00, 1.00)
    ),
    ncol = 2, nrow = 5
  )
  if (!identical(result, expected)) {
    print("ERROR")
    print("result:")
    print(result)
    print("expected:")
    print(expected)
  }
  expect_equal(identical(result, expected), TRUE)
})




test_that("stretch_nltt_matrix #3", {
  # Fill in the timepoints:
  #
  # t   N
  # 0.0 0.2
  # 0.4 0.5
  # 1.0 1.0
  #
  # becomes
  #
  # t   N
  # 0.0 0.2
  # 0.1 0.2
  # 0.2 0.2
  # 0.3 0.2
  # 0.4 0.5
  # 0.5 0.5
  # 0.6 0.5
  # 0.7 0.5
  # 0.8 0.5
  # 0.9 0.5
  # 1.0 1.0

  test <- matrix(c(c(0.0, 0.4, 1.0), c(0.2, 0.5, 1.0)), ncol = 2, nrow = 3)
  colnames(test) <- c("t", "N")
  expected <- matrix(
    c(
      seq(0.0, 1.0, 0.1),
      rep(0.2, times = 4), rep(0.5, times = 6), 1.0),
      ncol = 2, nrow = 11
  )
  result <- stretch_nltt_matrix(m = test, dt = 0.1, step_type = "lower")
  if (!identical(result, expected)) {
    print("ERROR")
    print("result:")
    print(result)
    print("expected:")
    print(expected)
  }
  expect_equal(identical(result, expected), TRUE)
})






test_that("get_average_nltt_matrix: stop on incorrect input", {
  test <- matrix(c(c(0.0, 0.4, 1.0), c(0.2, 0.5, 1.0)), ncol = 2, nrow = 3)
  colnames(test) <- c("t", "N")

  # must supply a matrx
  expect_error(
    stretch_nltt_matrix(m = list(), dt = 0.1, step_type = "upper")
  )

  # must supply at matrix with two columns
  expect_error(
    stretch_nltt_matrix(
      m = matrix(rep(0.0, times = 9), ncol = 3, nrow = 3),
      dt = 0.1,
      step_type = "upper"
    )
  )

  # step_type should be either 'upper' or 'lower'
  expect_silent(
    stretch_nltt_matrix(m = test, dt = 0.1, step_type = "upper")
  )
  expect_silent(
    stretch_nltt_matrix(m = test, dt = 0.1, step_type = "lower")
  )
  expect_error(
    stretch_nltt_matrix(m = test, dt = 0.1, step_type = "nonsense")
  )

})
