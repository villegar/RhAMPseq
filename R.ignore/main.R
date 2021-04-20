# Create Hex R-like logo
# hex_logo(dpi = 500, h_fill = "#4F0924", output = "figures/logo.png")
`%>%` <- magrittr::`%>%`

################################################################################
################################ Load raw data #################################
################################################################################
raw <- RhAMPseq::load_data(#fasta = "data/HaplotypeAllele.fasta",
                           hap_geno = "data/hap_genotype",
                           count_mat = "data/hap_genotype_matrix")

################################################################################
################################## Haplotypes ##################################
################################################################################
## Load groups by function
diversity <- readr::read_csv("data/diversity.csv", col_names = "name")
known <- readr::read_csv("data/known.csv", col_names = "name")
mapping <- readr::read_csv("data/mapping.csv", col_names = paste0("X", 1:7)) %>%
  .[, 1] %>%
  magrittr::set_names("name")
misc <- readr::read_csv("data/misc.csv", col_names = "name")
unknown <- readr::read_csv("data/unknown.csv", col_names = "name")

# 0. Create DB of known samples.
# 1. Find duplicates within replicates
# 2. Find unique fingerprint
# 3. Match from the known pool (across all of them)
hap_geno_mapping_only <- readr::read_tsv("data/hap_genotype")
hap_geno_all <- readr::read_tsv("data/old/hap_genotype")

## Clean haplotypes
hap2 <- hap_geno_mapping_only %>%
  dplyr::select(-Locus, -Haplotypes) %>%
  purrr::map_df(~purrr::map_chr(.x, RhAMPseq::cln_haplo, read_length = 5)) %>%
  dplyr::mutate(Locus = hap_geno_mapping_only$Locus, 
                Haplotypes = hap_geno_mapping_only$Haplotypes,
                .before = 1)

hap3 <- hap_geno_all %>%
  dplyr::select(-Locus, -Haplotypes) %>%
  purrr::map_df(~purrr::map_chr(.x, RhAMPseq::cln_haplo, read_length = 5)) %>%
  dplyr::mutate(Locus = hap_geno_all$Locus, 
                Haplotypes = hap_geno_all$Haplotypes,
                .before = 1)

# Reshape the tibbles, to long format
hap_known <- hap3 %>%
  dplyr::select(Locus, dplyr::any_of(known$name)) %>%
  tidyr::pivot_longer(c(dplyr::everything(), -Locus)) %>%
  dplyr::arrange(name) %>%
  dplyr::mutate(name = name %>% 
                  stringr::str_extract("[^__]*")) %>%
  dplyr::mutate(fp = FALSE)


hap_mapping <- hap2 %>%
  dplyr::select(Locus, dplyr::any_of(mapping$name)) %>%
  tidyr::pivot_longer(c(dplyr::everything(), -Locus)) %>%
  dplyr::arrange(name) %>%
  dplyr::mutate(name = name %>% 
                  stringr::str_extract("[^__]*")) %>%
  dplyr::mutate(fp = FALSE)

### Find "fingerprints"
#### Create reference: the parents, 588160 and 588271
reference <- hap_known %>%
  dplyr::filter(name %>% 
                  stringr::str_detect("^588160") |
                name %>% 
                  stringr::str_detect("^588271"))

aux <- c("588160") %>%
  purrr::map(function(r) {
    tmp <- reference %>%
      dplyr::filter(name %in% r) %>%
      dplyr::filter(!is.na(value))
    df <- unique(tmp$Locus) %>%
      purrr::map_df(function(l) {
        # print(l)
        tmp2 <- tmp %>%
          dplyr::filter(Locus %in% l) %>% 
          dplyr::mutate(spl = value %>%
                          purrr::map(split_haplo)) %>%
          tidyr::unnest(spl) %>%
          dplyr::rowwise() %>%
          dplyr::mutate(C = max(A, B),
                        A = min(A, B),
                        B = C) %>%
          dplyr::ungroup() %>%
          dplyr::count(Locus, name, A, B) %>%
          dplyr::mutate(freq = n / sum(n)) 
        # if (nrow(tmp2) == 2) {
        #   aux <- tmp2$A[2]
        #   tmp2$A <- tmp2$B[2]
        #   tmp2$B <- aux
        #   tmp3 <- tmp2 %>%
        #     dplyr::distinct(Locus, name, A, B)
        #   
        #   # tmp3 %>%
        #   #   dplyr::rowwise() %>%
        #   #   dplyr::select(A, B) %>%
        #   #   purrr::pmap(function(A, B) {
        #   #     tmp2 %>%
        #   #       dplyr::filter(A %in% A, B %in% B) %>%
        #   #       dplyr::select(n, freq) %>%
        #   #       colSums() %>%
        #   #       t() %>%
        #   #       tibble::as_tibble()
        #   #   })
        #     
        #   if (nrow(tmp3) == 1) {
        #     tmp2 <- tmp3 %>%
        #       dplyr::mutate(st = tmp2 %>%
        #                       dplyr::select(n, freq) %>%
        #                       colSums() %>%
        #                       t() %>%
        #                       tibble::as_tibble() %>%
        #                       list()) %>%
        #       tidyr::unnest(st)
        #   }
        #   #   tmp3 <- tmp3 %>%
        #   #   dplyr::mutate(n = su)
        # }
        tmp2
      })
    # l <- unique(tmp$Locus)
    # tmp2 <- tmp %>%
    #   dplyr::filter(Locus %in% l[1])
    # tmp2 
  })

hap_mapping$fp <-  unique(hap_mapping$name) %>%
  purrr::map(function(i) {
    s1 <- hap_mapping %>%
      dplyr::filter(name == i)
    s2 <- reference %>%
      dplyr::filter(name != i)
    !(s1$value %in% s2$value)
  }) %>%
  purrr::flatten_lgl()

hap_mapping_fp_mask <- hap_mapping %>%
  dplyr::select(-value) %>%
  tidyr::pivot_wider(Locus, names_from = name, values_from = fp)

hap_mapping_fp <- hap_mapping %>%
  dplyr::select(-fp) %>%
  tidyr::pivot_wider(Locus, names_from = name, values_from = value)

fp <- hap_mapping %>%
  dplyr::mutate(value = ifelse(fp, value, ifelse(is.na(value), value, "."))) %>%
  dplyr::select(-fp) %>%
  tidyr::pivot_wider(Locus, names_from = name, values_from = value)
readr::write_csv(fp, "~/Downloads/mapping-group-fingerprint.csv")

# This one is empty
hap_unknown <- hap2 %>%
  dplyr::select(Locus, dplyr::any_of(unknown$name)) %>%
  tidyr::pivot_longer(c(dplyr::everything(), -Locus)) %>%
  dplyr::arrange(name) %>%
  dplyr::mutate(fp = FALSE)

hap_mapping2 <- hap_mapping %>% 
  tidyr::pivot_longer(c(dplyr::everything(), -Locus)) %>%
  dplyr::arrange(name) %>%
  dplyr::mutate(fp = FALSE)

## Load mapping family genotypes
mapping_family <- RhAMPseq::read_excel_col("data/Rhampseq_populations.xlsx", "A")
mapping_family_short <- unlist(lapply(mapping_family, 
                                      function(x) gsub("__.*", "", x)))
mapping_family_short <- unique(mapping_family_short)

## Extract genotypes from the hap_genotype file
hap_geno_names <- names(raw$hap_geno)
hap_geno_names_blank <- hap_geno_names[grepl("*blank*", hap_geno_names)]
## Drop "blank" genotypes
# hap_geno_names <- hap_geno_names[!grepl("*blank*", hap_geno_names)]
## Trim 'X' from genotypes [added when importing the data]
hap_geno_names <- gsub("X*", "", hap_geno_names)
## Further trimming of genotypes: 160_271_001__vDNAlon499A03_A01 -> 160_271_001
hap_geno_names <- unlist(lapply(hap_geno_names, 
                                function(x) gsub("__.*", "", x)))

## Lookup for mapping_family genotypes
hap_geno_names_idx <- hap_geno_names %in% mapping_family_short
## Extract matching genotypes
hap_geno_names_sub <- hap_geno_names[hap_geno_names_idx]
hap_geno_sub <- raw$hap_geno[hap_geno_names_idx]

## Extract locus data
locus <- raw$hap_geno$Locus

# Start parallel backend
tictoc::tic("Parallel") # ~12.719 sec (2 CPUs); ~7.923 sec (4 CPUs)
cores <- parallel::detectCores() - 1
cpus <- ifelse(cores > 2, 2, 1)
cl <- parallel::makeCluster(cpus, setup_strategy = "sequential")
doParallel::registerDoParallel(cl)

# Load binary operator for backend
`%dopar%` <- foreach::`%dopar%`
map <- foreach::foreach(i = 1:length(hap_geno_sub), 
                        .combine = cbind) %dopar% {
                          unlist(lapply(hap_geno_sub[, i], get_geno))
                        }
parallel::stopCluster(cl) # Stop cluster
tictoc::toc()
colnames(map) <- names(hap_geno_sub)
rownames(map) <- locus
# map <- cbind(Locus = locus, map)
map_t <- t(map)
colnames(map_t) <- col(map)
rownames(map_t) <- gsub("X*", "", rownames(map_t))
write.csv(map_t, "data/pseudomap.csv")
write.csv(map_t[1:5,1:5], "data/pseudomap-sample.csv")

# tictoc::tic("Serial") # ~21.387 sec
# map_s <- data.frame(Locus = locus)
# for(i in 1:length(hap_geno_sub)) {
#   map_s <- cbind(map_s, unlist(lapply(hap_geno_sub[, i], get_geno)))
#   colnames(map_s)[i + 1] <- names(hap_geno_sub)[i] # i + 1 column
# }
# tictoc::toc()

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