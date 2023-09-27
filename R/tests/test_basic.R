helpers_dir <- "/n/data1/dfci/pedonc/kadoch/ben/soft_ben/kadoch_pipeline_tools/R/rna_seq"
r_dir <- function(x){source(file.path(helpers_dir, x))}

library("testthat")
source(rdir("dge_helpers.R"))
source(rdir('additive.R'))


adjpval_threshold <- 0.05
logfc_threshold <- 1

raw_counts <- read.csv("tests/dummy_reads.csv")
names(raw_counts)[1] <- 'gene_names'

cpms <- log2(raw_counts[,2:ncol(raw_counts)])
cpms <- cbind(raw_counts[1], cpms)

groups <- c(rep('control', 3), rep('t1', 3), rep('t2', 3), rep('t3', 3))

# Test 1
test_that("DGE directionality is correct", {
  
  de <- get_de_between_conditions(cpms, "control", "t1", groups, 
                                  col_offset=1) 
  
  expect_equal(de$logFC[de$gene_name == 'DEEZ'], -4)
  
  de <- get_de_between_conditions(cpms, "t1", "control", groups, 
                                  col_offset=1)
  
  expect_equal(de$logFC[de$gene_name == 'DEEZ'], 4)
  
  
  de <- get_de_between_conditions(cpms, "control", "t2", groups, 
                                  col_offset=1)
  expect_equal(get_de_gene_names(de, 0.5, 0.5), character())
  
  de <- get_de_between_conditions(cpms, "t2", "t1", groups, 
                                  col_offset=1) 
  
  expect_equal(de$logFC[de$gene_name == 'DEEZ'], -4)
  
  de <- get_de_between_conditions(cpms, "t1", "t2", groups, 
                                  col_offset=1)
  expect_equal(de$logFC[de$gene_name == 'DEEZ'], 4)
  
})

test_that('DE genes called correctly',{
  # Test absolute DE gene regulation
  
  de_binary <- get_de_binary(cpms, 'control', 't1', groups, adjpval_threshold,
                             logfc_threshold)
  
  de <- get_de_between_conditions(cpms, "control", "t1", groups, 
                                  col_offset=1) 
  
  de_gene_names <- get_de_gene_names(de, adjpval_threshold, logfc_threshold)
  
  expect_equal(de_gene_names, names(de_binary)[de_binary], 
               c("DEEZ",  "JEFF",  "DRAKE"))
  
  # Test positive regulation
  de_binary <- get_de_binary(cpms, 'control', 't1', groups, adjpval_threshold,
                             logfc_threshold, thr_set='increase')
  
  de <- get_de_between_conditions(cpms, "control", "t1", groups, 
                                  col_offset=1) 
  
  de_gene_names <- get_de_gene_names(de, adjpval_threshold, logfc_threshold, thr_set='increase')
  
  expect_equal(de_gene_names, names(de_binary)[de_binary], 
               c("JEFF",  "DRAKE"))
  
  # Test decreased expression
  de_binary <- get_de_binary(cpms, 'control', 't1', groups, adjpval_threshold,
                             logfc_threshold, thr_set='decrease')
  
  de <- get_de_between_conditions(cpms, "control", "t1", groups, 
                                  col_offset=1) 
  
  de_gene_names <- get_de_gene_names(de, adjpval_threshold, logfc_threshold, thr_set='decrease')
  
  expect_equal(de_gene_names, names(de_binary)[de_binary], 
               c("DEEZ"))
  
})

undebug(get_synergy_matrix)
undebug(get_mean_expression)

test_that("Synergy matrix is created properly.",{
  
  raw_counts <- read.csv("tests/additive_reads.csv")
  names(raw_counts)[1] <- 'gene_names'
  
  cpms <- log2(raw_counts[,2:ncol(raw_counts)])
  cpms <- cbind(raw_counts[1], cpms)
  
  groups <- c(rep('control', 3), rep('t1', 3), rep('t2', 3), rep('t3', 3))
  
  synergy_matrix <- get_synergy_matrix(cpms, "control", 't1', 't2', "t3",
                                       groups, num_desc_columns=1)
  
  expect_equal(synergy_matrix['SOUP', 'S'], -0.25)
  expect_equal(synergy_matrix['DRAKE', 'S'], 0.1)
  expect_equal(synergy_matrix['JEFF', 'S'], 4)
  expect_equal(synergy_matrix['DEEZ', 'S'], 0)
  
})


