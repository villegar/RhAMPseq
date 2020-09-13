hex_logo <- function(subplot = "figures/slope.png",
                     dpi = 800,
                     h_color = "#000000",
                     h_fill = "#363b74",
                     output = "figures/metapipe.png",
                     p_color = "#eeeeee",
                     url = "https://github.com/villegar/RhAMPseq") {
  hexSticker::sticker(subplot = subplot, package = "RhAMPseq", 
                      h_color = h_color,  h_fill = h_fill,
                      dpi = dpi,
                      s_x = 1.0, s_y = .85, s_width = .5,
                      p_x = 1.0, p_y = 1.52, p_size = 6, p_color = p_color,
                      url = url, 
                      u_angle = 30, u_color = p_color, u_size = 1.35,
                      filename = output)
}

load_data <- function(fasta = NULL, hap_geno = NULL, count_mat = NULL) {
  fasta_data <- NULL
  hap_data <- NULL
  count_mat_data <- NULL
  
  if(!is.null(fasta)) {
    fasta_data <- seqinr::read.fasta(fasta)
  }
  
  if(!is.null(hap_geno)) {
    hap_data <- read.table(hap_geno, header = TRUE) 
  }
  
  if(!is.null(count_mat)) {
    count_mat_data <- read.table(count_mat, header = TRUE) 
  }
  
  return(list(fasta = fasta_data, hap_geno = hap_data, count_mat = count_mat_data))
}

get_seq <- function(raw, upper = TRUE) {
  seq <- paste0(unname(unlist(raw)), collapse = "")
  return(ifelse(upper, toupper(seq), seq))
}


#' Read columns from MS Excel file
#'
#' @param filename filename, absolute or relative paths are valid
#' @param columns vector with the columns names
#'
#' @return table with the requested columns
# @export
#'
#' @examples
#' read_excel_col("data/Rhampseq_populations.xlsx", "A")
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
# @export
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