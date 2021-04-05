#' Create hexagonal logo for the package
#'
#' @param subplot image to use as the main logo
#' @param dpi plot resolution (dots-per-inch)
#' @param h_color color for hexagon border
#' @param h_fill color to fill hexagon
#' @param output output file (hexagonal logo)
#' @param package title for logo (package name)
#' @param p_color color for package name
#' @param url URL for package repository or website
#' @param u_size text size for URL
#'
#' @return hexagonal logo
#' @export
#'
#' @examples
#' hex_logo()
#' \dontrun{
#' hex_logo("inst/images/slope.png", output = "inst/images/logo.png")
#' }
hex_logo <- function(subplot = system.file("images/slope.png", package = "RhAMPseq"),
                     dpi = 600,
                     h_color = "#000000",
                     h_fill = "#363B74",
                     output = system.file("images/logo.png", package = "RhAMPseq"),
                     p_color = "#EEEEEE",
                     url = "https://github.com/villegar/RhAMPseq",
                     u_size = 1.35) {
  hexSticker::sticker(subplot = subplot, package = "RhAMPseq", 
                      h_color = h_color,  h_fill = h_fill,
                      dpi = dpi,
                      s_x = 1.0, s_y = .85, s_width = .5,
                      p_x = 1.0, p_y = 1.52, p_size = 6, p_color = p_color,
                      url = url, 
                      u_angle = 30, u_color = p_color, u_size = u_size,
                      filename = output)
}

#' Load rhAmpSeq raw data in different formats
#'
#' @param fasta FASTA data 
#' @param hap_geno genotypes data
#' @param count_mat count matrix
#'
#' @return list of each loaded file
#' @export
#'
#' @examples
#' \dontrun{
#' load_data(fasta = "data/HaplotypeAllele.fasta",
#'           hap_geno = "data/hap_genotype",
#'           count_mat = "data/readCountMatrixFile")
#' }
load_data <- function(fasta = NULL, hap_geno = NULL, count_mat = NULL) {
  fasta_data <- NULL
  hap_data <- NULL
  count_mat_data <- NULL
  
  if(!is.null(fasta)) {
    fasta_data <- seqinr::read.fasta(fasta)
  }
  
  if(!is.null(hap_geno)) {
    hap_data <- readr::read_tsv(hap_geno, col_names = TRUE) 
  }
  
  if(!is.null(count_mat)) {
    count_mat_data <- readr::read_tsv(count_mat, col_names = TRUE) 
  }
  
  return(list(fasta = fasta_data, hap_geno = hap_data, count_mat = count_mat_data))
}

#' Combine list of bases
#'
#' @param raw list of bases
#' @param upper flag on whether the sequence should be returned in all caps
#'
#' @return string of bases
#' @export
#'
#' @examples
#' set.seed(123)
#' # Generate random list of DNA bases
#' seq_list <- list(c("A", "C", "G", "T")[sample(1:4, 100, TRUE)])
#' get_seq(seq_list) 
#' get_seq(seq_list, FALSE) # Return lower case sequence
get_seq <- function(raw, upper = TRUE) {
  seq <- paste0(unname(unlist(raw)), collapse = "")
  return(ifelse(upper, toupper(seq), tolower(seq)))
}


#' Read columns from MS Excel file
#'
#' @param filename filename, absolute or relative paths are valid
#' @param columns vector with the columns names
#'
#' @return table with the requested columns
#' @export
#'
#' @examples
#' \dontrun{
#' read_excel_col("data/Rhampseq_populations.xlsx", "A")
#' }
read_excel_col <- function(filename, columns) {
  return(readxl::read_excel(filename, range = readxl::cell_cols(columns)))
}

#' Get genotypes from haplotypes
#'
#' @param hap haplotypes data with format: X/Y:n,m
#' @param genotypes genotypes to assign based on homozygosity, 
#' c(Homozygous, Heterozygous, NULL); NULL haplotypes [./.:0]
#'
#' @return genotype
#' @export
#'
#' @examples
#' get_geno("./.:0")
#' get_geno("4/2:34,24")
#' get_geno("2/2:71")
get_geno <- function(hap, genotypes = c("NN", "NP", NA)) {
  if(hap == "./.:0")
    return(genotypes[3])
  # Split haplotypes and read count
  hap_arr <- unlist(strsplit(hap, ":"))
  # Extract haplotypes
  hap_names <- unlist(strsplit(hap_arr[1], "/"))
  ## Compare if variance = 0 (Homozygous)
  return(ifelse(var(hap_names) == 0, genotypes[1], genotypes[2]))
  
  # Extract read count
  # unlist(strsplit(hap_arr[2], ","))
}

#' Clean haplotypes
#' 
#' Search for haplotypes below a threshold for read length and set them to
#' \code{missing}.
#'
#' @param hap haplotype data with the format: \code{X/Y:n,m}.
#' @param read_length minimum read length to keep haplotype.
#' @param missing default value for haplotypes below the \code{read_length}.
#' 
#' @return haplotype
#' @export
#'
#' @examples
#' cln_haplo("./.:0")
#' cln_haplo("4/2:34,24")
#' cln_haplo("2/2:71")
#' cln_haplo("4/2:2,1")
#' cln_haplo("2/2:4")
cln_haplo <- function(hap, read_length = 5, missing = NA) {
  hap <- trimws(hap)
  if(hap == "./.:0")
    return(missing)
  # Split haplotypes and read count
  hap_arr <- unlist(strsplit(hap, ":"))
  # Extract haplotypes
  hap_names <- unlist(strsplit(hap_arr[1], "/"))
  # ## Compare if variance = 0 (Homozygous)
  # return(ifelse(var(hap_names) == 0, genotypes[1], genotypes[2]))
  
  # Extract read count
  tryCatch({
    rl <- unlist(strsplit(hap_arr[2], ",")) %>%
      purrr::map_dbl(as.numeric) %>%
      sum(na.rm = TRUE)
    if (rl < read_length)
      return(missing)
  }, error = function(e) {
    stop(e)
  })
  return(hap)
}
