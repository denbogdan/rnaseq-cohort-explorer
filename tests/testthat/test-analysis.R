source("../../R/data.R"); source("../../R/analysis.R")

set.seed(1)
mat <- matrix(rexp(50 * 40, 1), nrow = 50, dimnames = list(paste0("G", 1:50), paste0("P", 1:40)))
mat["G2", ] <- mat["G1", ] * 2 + rnorm(40, sd = 0.01)
surv <- data.frame(patient_id = colnames(mat), os_months = rexp(40, 0.05), death = rbinom(40, 1, 0.6))

test_that("correlate_genes finds a monotone relationship", {
  r <- correlate_genes(mat, "G1", "G2")
  expect_gt(r$rho, 0.95)
  expect_lt(r$p, 0.001)
})

test_that("unknown genes raise a readable error", {
  expect_error(correlate_genes(mat, "G1", "NOPE"), "not found")
})

test_that("rank_by_correlation is decreasing and drops the query gene", {
  rk <- rank_by_correlation(mat, "G1", rownames(mat))
  expect_false("G1" %in% names(rk))
  expect_equal(names(rk)[1], "G2")
  expect_false(is.unsorted(rev(rk)))
})

test_that("median stratification yields two groups", {
  s <- stratify_survival(mat, surv, "G1", maxstat = FALSE)
  expect_setequal(unique(s$data$group), c("high", "low"))
})

test_that("top_expressed returns the requested number of genes", {
  expect_length(top_expressed(mat, 10), 10)
})
