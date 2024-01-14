# test_input_processing.R
library(testthat)


source("src/input_processing.R")

# Test that read_htseq_files function reads files correctly
test_that("read_htseq_files reads files correctly", {
  # Create temporary directory
  temp_dir <- tempdir()
  
  # Create test data frames
  df1 <- data.frame(V1 = c("A", "B", "C"), V2 = c(10, 20, 30))
  df2 <- data.frame(V1 = c("D", "E", "F"), V2 = c(40, 50, 60))
  
  # Write test data frames to temporary files
  write.table(df1, file.path(temp_dir, "file1.txt"), sep = "\t", row.names = FALSE)
  write.table(df2, file.path(temp_dir, "file2.txt"), sep = "\t", row.names = FALSE)
  
  # Call function with path to temporary directory
  result <- read_htseq_files(temp_dir)
  
  # Check results
  expect_equal(length(result), 2)
  expect_equal(names(result), c("file1", "file2"))
  expect_equal(result$file1, df1)
  expect_equal(result$file2, df2)
})

# Test that process_ht_seq_table function processes table correctly
test_that("process_ht_seq_table processes table correctly", {
  input <- data.frame(V1 = c("A", "B", "C"), V2 = c(10, 20, 30))
  result <- process_ht_seq_table(input)
  expect_equal(result, data.frame(V1 = c("Top5", "Not Top5"), SUM = c(30, 30)))
})

# Test that create_sample_table function creates table correctly
test_that("create_sample_table creates table correctly", {
  # Mock rbindlist and pivot_wider functions to simulate table creation
  with_mock(
    rbindlist = function(idcol) data.frame(V1 = c("A", "B", "C"), Sample = c("file1", "file2", "file3"), V2 = c(10, 20, 30)),
    pivot_wider = function(names_from, values_from) data.frame(V1 = c("A", "B", "C"), file1 = c(10, 20, 30), file2 = c(10, 20, 30), file3 = c(10, 20, 30)),
    {
      result <- create_sample_table(list("file1.txt", "file2.txt", "file3.txt"))
      expect_equal(result, data.frame(V1 = c("A", "B", "C"), file1 = c(10, 20, 30), file2 = c(10, 20, 30), file3 = c(10, 20, 30)))
    }
  )
})

# Test that combine_tpms function combines TPM files correctly
test_that("combine_tpms combines TPM files correctly", {
  # Mock list.files and fread functions to simulate file reading
  with_mock(
    list.files = function(pattern, recursive, full.names) c("file1.txt", "file2.txt"),
    fread = function(file) data.frame(V1 = c("A", "B", "C"), V2 = c(10, 20, 30)),
    {
      result <- combine_tpms("/path/to/files")
      expect_equal(result, data.frame(GeneID = c("A", "B", "C"), file1 = c(10, 20, 30), file2 = c(10, 20, 30)))
    }
  )
})

# Test that calculate_z_scores function calculates z-scores correctly
test_that("calculate_z_scores calculates z-scores correctly", {
  input <- data.frame(GeneID = c("A", "B", "C"), file1 = c(10, 20, 30), file2 = c(10, 20, 30))
  result <- calculate_z_scores(input)
  expect_equal(result, matrix(c(-1, 0, 1, -1, 0, 1), nrow = 3, ncol = 2))
})