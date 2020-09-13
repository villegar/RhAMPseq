test_that("hexagonal logo works", {
  hex_logo(output = "hex_logo.png")
  expect_true(file.exists("hex_logo.png"))
  expect_false(dir.exists("hex_logo.png"))
  expect_gt(file.size("hex_logo.png"), 0)
  file.remove("hex_logo.png")
  expect_false(file.exists("hex_logo.png"))
})

test_that("load data works", {
  # Create test FASTA file
  set.seed(123)
  seq001 <- paste0(c("A", "C", "G", "T")[sample(1:4, 100, TRUE)], collapse = "")
  seq002 <- paste0(c("A", "C", "G", "T")[sample(1:4, 100, TRUE)], collapse = "")
  seqinr::write.fasta(sequences = list(seq001, seq002), 
                      names = list("001", "002"), 
                      file.out = "test.fasta")
  # Create test genotypes file
  hap_data <- data.frame(A = 1:10,
                         B = NA,
                         C = letters[1:10])
  write.table(hap_data, "test.hap")
  # Create test count matrix
  count_mat <- data.frame(A = 1:10,
                          B = NA,
                          C = letters[1:10])
  write.table(count_mat, "test.count_mat")
  
  # Check for test files
  filenames <- c("test.fasta",
                 "test.hap",
                 "test.count_mat")
  for (f in filenames) {
    expect_true(file.exists(f))
    expect_false(dir.exists(f))
    expect_gt(file.size(f), 0)
  }
  
  raw <- load_data(fasta = "test.fasta",
                   hap_geno = "test.hap",
                   count_mat = "test.count_mat")
  
  expect_equal(length(raw), 3)
  expect_equal(names(raw), c("fasta", "hap_geno", "count_mat"))
  
  # Delete files
  for (f in filenames) {
    file.remove(f)
    expect_false(file.exists(f))
  }
})