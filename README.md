
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RhAMPseq: Analysis of rhAmpSeq data <img src="https://raw.githubusercontent.com/villegar/RhAMPseq/master/inst/images/logo.png" alt="logo" align="right" height=200px/>

<!-- Analysis of rhAmpSeq data. -->

<!-- badges: start -->

[![R build
status](https://github.com/villegar/RhAMPseq/workflows/R-CMD-check/badge.svg)](https://github.com/villegar/RhAMPseq/actions)
[![](https://img.shields.io/badge/devel%20version-0.0.1-blue.svg)](https://github.com/villegar/MetaPipe)
[![](https://codecov.io/gh/villegar/RhAMPseq/branch/master/graph/badge.svg)](https://codecov.io/gh/villegar/RhAMPseq)
<!-- badges: end -->

## Installation

<!-- You can install the released version of RhAMPseq from [CRAN](https://CRAN.R-project.org) with: -->

<!-- ``` r -->

<!-- install.packages("RhAMPseq") -->

<!-- ``` -->

<!-- And the development version from [GitHub](https://github.com/) with: -->

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages(c("hexSticker", "remotes")
remotes::install_github("villegar/RhAMPseq")
```

## Data description

**Author:**
    J.L.

### Files

  - <span style="color: #03045e; background-color: #E5F1FF; padding: 2px 5px;">hap\_genotype</span>:
    haplotype genotyping data for all samples. This is in a number code.
    The number code corresponds to a specific sequenced haplotype at
    each
    locus.
  - <span style="color: #03045e; background-color: #E5F1FF; padding: 2px 5px;">HaplotypeAllele.fasta</span>:
    each marker is a conserved region of the genome, so the haplotype
    allele may differ only by a single SNP, or it could be many SNPS or
    indels.  
  - <span style="color: #03045e; background-color: #E5F1FF; padding: 2px 5px;">readCountMatrixFile</span>:
    number of read for each amplicon and each sample. This data is also
    provided in the
    <span style="color: #03045e; background-color: #E5F1FF; padding: 2px 5px;">hap\_genotype</span>
    file.

### Example

| Locus              | Haplotypes                                                                        | blank\_\_vDNAfen551A01\_A01 | 160\_271\_303\_\_vDNAfen551A01\_A02 | 160\_271\_311\_\_vDNAfen551A01\_A03 |
| :----------------- | :-------------------------------------------------------------------------------- | :-------------------------- | :---------------------------------- | :---------------------------------- |
| rh\_chr1\_10143990 | 2(0.37);4(0.16);9(0.10);5(0.09);1(0.07);3(0.07);8(0.05);13(0.04);6(0.03);7(0.03); | ./.:0                       | 4/2:34,24                           | 2/2:71                              |

where:

  - `Locus` is the marker.

  - `Haplotypes` lists the frequency of haplotypes at this `Locus` in
    the *ENTIRE* dataset, not just our samples. You can mostly ignore
    this.

  - Next column starts with `blank`, this is a blank control well in the
    plate. Here at this `Locus` it sequenced as `./.`, which means
    `NULL`. and `:0` which is the read depth.

  - Next column is sample `303` from our mapping family, its on the
    `fen551` plate. Here you see we detect `4/2` and `34, 24` for
    haplotype 4 with 34 reads, and 2 with 24 reads. This locus is
    *heterozygous*.

  - Next column is sample `311` from our mapping family. Here you see we
    only detect haplotype `2` at `71` reads. This is *homozygous*.

So, if we wanted to know the DNA sequence of this marker, we would look
up haplotypes 4 and 2 from the
<span style="color: #03045e; background-color: #E5F1FF; padding: 2px 5px;">HaplotypeAllele.fasta</span>
file.

One other thing to note that is very important. the marker name,
`rh_ch1_10143990` stands for rhampseq, chromosome 1, position 10143990.
This position information does not indicate the position of the
polymorphism detected at this marker, it is the position of the primer
sequence. So when we map these markers, they will be very near, but not
exactly at, these positions.

## Initial analysis

### Raw data

Load raw data from the [3 files](#files) previously described.

``` r
raw <- load_data(fasta = "data/HaplotypeAllele.fasta", 
                 hap_geno = "data/hap_genotype",
                 count_mat = "data/readCountMatrixFile")
```

### Haplotypes

``` r
# Load haplotypes for mapping family
mapping_family <- read_excel_col("data/Rhampseq_populations.xlsx", "A")
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
# map_t[1:5,1:5]
```

#### Pseudo-genetic map sample (`map_t`)

| X                                   | rhMAS\_5GT\_cons95 | rhMAS\_SDI\_p2\_AG11\_chr18\_30Mb | rhMAS\_SDI\_p3\_AGL11 | rhMAS\_SDI\_p4\_AGL11 | rhMAS\_5GT\_700 |
| :---------------------------------- | :----------------- | :-------------------------------- | :-------------------- | :-------------------- | :-------------- |
| 160\_271\_303\_\_vDNAfen551A01\_A02 | NP                 | NP                                | NP                    | NP                    | NN              |
| 160\_271\_311\_\_vDNAfen551A01\_A03 | NN                 | NN                                | NP                    | NP                    | NN              |
| 160\_271\_319\_\_vDNAfen551A01\_A04 | NP                 | NP                                | NN                    | NN                    | NN              |
| 160\_271\_327\_\_vDNAfen551A01\_A05 | NP                 | NP                                | NP                    | NP                    | NN              |
| 160\_271\_335\_\_vDNAfen551A01\_A06 | NP                 | NP                                | NP                    | NP                    | NN              |

### Other stuff

  - Extract marker names and drop their “repeat ID”, `#X`. e.g.
    `rhMAS_5GT_cons95#1` -\> `rhMAS_5GT_cons95`

<!-- end list -->

``` r
marker_names <- names(raw$fasta)
marker_names_no_repeats <- lapply(marker_names, function(x) gsub("#.*", "", x))

length(marker_names_no_repeats) # 984030
length(unique(marker_names_no_repeats)) # 2055
```

  - *Something that might not make sense …*

<!-- end list -->

``` r
sum_stats <- list() #data.frame(marker = NA, non_zero = NA)
for (loc in raw$hap_geno$Locus) {
  # repeats_idx <- which(marker_names_no_repeats == loc)
  # repeats <- raw$fasta[repeats_idx]
  # repeats_seq <- lapply(repeats, get_seq)
  
  count_mat_idx <- which(loc == names(raw$count_mat))
  count_mat_sub <- raw$count_mat[count_mat_idx]
  sum_stats[[loc]] <- list(idx = which(count_mat_sub > 0), count = length(unlist(count_mat_sub)))
  #sum_stats <- rbind(sum_stats, c(loc, length(which(count_mat_sub > 0))))
}
# sum_stats[1, ] <- NULL
```
