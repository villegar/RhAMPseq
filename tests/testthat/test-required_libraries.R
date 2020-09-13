test_that("check libraries function works", {
  libs <- c("doParallel",
            "foreach",
            "hexSticker",
            "knitr", 
            "parallel",
            "readxl", 
            "seqinr",
            "tictoc")
  if(!require(l, character.only = TRUE)) {
    remove.packages("tictoc")
  }
  check_libs()
    
  for(l in libs) {
    expect_true(require(l, character.only = TRUE))
  }
})
