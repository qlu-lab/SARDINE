sigma_tsv <- system.file("extdata", "fixtures", "inputs", "sigma_baseline_synthetic.tsv", package = "sardine")

test_that("read_dnam_covariance accepts a file, matrix, data.frame and data.table", {
  from_file <- read_dnam_covariance(sigma_tsv, verbose = FALSE)
  expect_true(is.matrix(from_file))
  expect_identical(rownames(from_file), colnames(from_file))

  expect_identical(read_dnam_covariance(from_file, verbose = FALSE), from_file)
  expect_identical(read_dnam_covariance(as.data.frame(from_file), verbose = FALSE), from_file)

  with_id_col <- data.table::data.table(CpGmarker = rownames(from_file), from_file)
  expect_identical(read_dnam_covariance(with_id_col, verbose = FALSE), from_file)

  rds <- tempfile(fileext = ".rds")
  saveRDS(from_file, rds)
  expect_identical(read_dnam_covariance(rds, verbose = FALSE), from_file)
})

test_that("read_dnam_covariance rejects mismatched names, text columns and other types", {
  m <- matrix(c(1, 0.5, 0.5, 1), 2, dimnames = list(c("a", "b"), c("b", "a")))
  expect_error(read_dnam_covariance(m, verbose = FALSE), "identical and in the same order")

  bad <- data.table::data.table(id = c("a", "b"), a = c("x", "y"), b = c(1, 2))
  expect_error(read_dnam_covariance(bad, verbose = FALSE), "non-numeric covariance columns: a")

  expect_error(read_dnam_covariance(1:3, verbose = FALSE), "must be a delimited file")
  expect_error(read_dnam_covariance("does_not_exist.tsv", verbose = FALSE), "File not found")
})
