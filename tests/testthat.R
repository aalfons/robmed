if(requireNamespace("testthat", quietly = TRUE)) {
  library("robmed", quietly = TRUE)
  testthat::test_check("robmed")
}
