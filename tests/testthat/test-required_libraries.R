test_that("check libraries function works", {
  check_libs()
  libs <- c("doParallel",
            "foreach",
            "hexSticker",
            "knitr", 
            "parallel",
            "readxl", 
            "seqinr",
            "tictoc") 
  for(l in libs) {
    expect_true(require(l, character.only = TRUE))
  }
})
