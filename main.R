# Create Hex R-like logo
# hex_logo(dpi = 500, h_fill = "#4F0924", output = "figures/logo.png")

# Load raw data
raw <- load_data(fasta = "data/HaplotypeAllele.fasta", 
                 hap_geno = "data/hap_genotype",
                 count_mat = "data/readCountMatrixFile")

################################################################################
################################## Haplotypes ##################################
################################################################################
## Load mapping family genotypes
mapping_family <- read_excel_col("data/Rhampseq_populations.xlsx", "A")
mapping_family_short <- unlist(lapply(mapping_family, function(x) gsub("__.*", "", x)))

## Extract genotypes from the hap_genotype file
hap_geno_names <- names(raw$hap_geno)
hap_geno_names_blank <- hap_geno_names[grepl("*blank*", hap_geno_names)]
## Drop "blank" genotypes
# hap_geno_names <- hap_geno_names[!grepl("*blank*", hap_geno_names)]
## Trim 'X' from genotypes [added when importing the data]
hap_geno_names <- gsub("X*", "", hap_geno_names)
## Further trimming of genotypes: 160_271_001__vDNAlon499A03_A01 -> 160_271_001
hap_geno_names <- unlist(lapply(hap_geno_names, function(x) gsub("__.*", "", x)))

## Lookup for mapping_family genotypes
hap_geno_names_idx <- hap_geno_names %in% mapping_family_short
## Extract matching genotypes
hap_geno_names_sub <- hap_geno_names[hap_geno_names_idx]
hap_geno_sub <- raw$hap_geno[hap_geno_names_idx]

################################################################################
############################## Reads Count Matrix ##############################
################################################################################
sum_stats <- data.frame(marker = NA, non_zero = NA)
for (loc in raw$hap_geno$Locus) {
  count_mat_idx <- which(loc == names(raw$count_mat))
  count_mat_sub <- raw$count_mat[count_mat_idx]
  sum_stats <- rbind(sum_stats, c(loc, length(which(count_mat_sub > 0))))
}
sum_stats <- sum_stats[-1, ]  # Drop first row

################################################################################
#################################### FASTA #####################################
################################################################################
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
}