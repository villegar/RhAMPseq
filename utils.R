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
  
  return(list(fasta_data, hap_data, count_mat_data))
}
