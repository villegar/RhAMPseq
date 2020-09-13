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

test_that("get sequence works", {
  set.seed(123)
  # Generate random list of DNA bases
  seq_list <- list(c("A", "C", "G", "T")[sample(1:4, 100, TRUE)])
  expect_equal(get_seq(seq_list), paste0(unlist(seq_list), collapse = ""))
  # Return lower case sequence
  expect_equal(get_seq(seq_list, FALSE), 
               tolower(paste0(unlist(seq_list), collapse = "")))
})

test_that("read Excel column works", {
  # Create test Excel data
  excel_data <- data.frame(A = 1:10,
                           B = NA,
                           C = letters[1:10])
  xlsx::write.xlsx(excel_data, "test.xls", row.names = FALSE)
  # Test that the Excel file was created
  expect_true(file.exists("test.xls"))
  expect_false(dir.exists("test.xls"))
  expect_gt(file.size("test.xls"), 0)
  
  col_data <- read_excel_col("test.xls", c("A"))
  expect_equal(nrow(col_data), 10)
  expect_equal(col_data$A, 1:10)
  
  # Delete test Excel file
  file.remove("test.xls")
  expect_false(file.exists("test.xls"))
})

test_that("get genotype works", {
  expect_true(is.na(get_geno("./.:0")))
  expect_equal(get_geno("4/2:34,24"), "NP")
  expect_equal(get_geno("2/2:71"), "NN")
  expect_error(get_geno(NULL))
})