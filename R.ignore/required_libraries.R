#' Check for required libraries and install if not currently available 
#'
#' @export
#'
#' @examples
#' \dontrun{
#' check_libs()
#' }
check_libs <- function() {
  libs <- c("doParallel",
            "foreach",
            "hexSticker",
            "knitr", 
            "magrittr",
            "parallel",
            "readxl", 
            "seqinr",
            "tictoc") 
  
  out <- lapply(libs, 
                function(x) {
                  if (!require(x, character.only = TRUE)) install.packages(x, dependencies = TRUE)
                })
}