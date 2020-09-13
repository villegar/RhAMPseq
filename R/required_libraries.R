libs <- c("doParallel",
          "foreach",
          "hexSticker",
          "knitr", 
          "parallel",
          "readxl", 
          "seqinr",
          "tictoc") 

out <- lapply(libs, 
       function(x) {
         if (!require(x, character.only = TRUE))
           install.packages(x, dependencies = TRUE)
       })