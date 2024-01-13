# Load the testthat package
library(testthat)

# print a message to the console
message("Running tests...")

# Load your script
source("src/dge.R")

# Start your test file
test_that("Test normalize function", {

  print("Running tests for normalize function...")

# Create a dummy DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = matrix(rpois(100, lambda = 10), 10),
                              colData = data.frame(condition = rep(c("A", "B"), each = 5)),
                              design = ~ condition)
  
  # Test the normalize function
  normalized_counts <- normalize(dds)
  
  # Check if the output is a matrix
  expect_is(normalized_counts, "matrix")
  
  # Check if the output has the correct dimensions
  expect_equal(dim(normalized_counts), dim(counts(dds)))
})

test_that("Test read_sample_files function", {

  #print a message to the console
  message("Running tests for read_sample_files function...")
  
  # Test the read_sample_files function with example files
  result <- read_sample_files("data/BKUS_SAMPLES/feature_table.csv", "data/BKUS_SAMPLES/metadata.tsv")
  
  # Check if the output is a list
  expect_is(result, "list")
  
  # Check if the list has the correct names
  expect_equal(names(result), c("feature_table", "metadata"))
})

# Add more tests for other functions in your script