# Create Hex R-like logo
# hex_logo(dpi = 500, h_fill = "#4F0924", output = "figures/logo.png")

# Load raw data
raw <- load_data(fasta = "data/HaplotypeAllele.fasta", 
                 hap_geno = "data/hap_genotype",
                 count_mat = "data/readCountMatrixFile")

marker_names <- names(raw$fasta)
# tictoc::tic("lapply")
marker_names_no_repeats <- lapply(marker_names, function(x) gsub("#.*", "", x))
# tictoc::toc()

length(marker_names_no_repeats)
length(unique(marker_names_no_repeats))

for (loc in raw$hap_geno$Locus) {
  repeats_idx <- which(marker_names_no_repeats == loc)
  repeats <- raw$fasta[repeats_idx]
  repeats_seq <- lapply(repeats, get_seq)
  
  count_mat_idx <- which(loc == names(raw$count_mat))
  count_mat_sub <- raw$count_mat[count_mat_idx]
}
