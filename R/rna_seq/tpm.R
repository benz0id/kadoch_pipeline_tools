source('/Users/btudorpr/Desktop/r_codebase/gene_mapping.R')


s <- function(x){as.character(x)}

tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

rpkm <- function(counts,len) {
  cpm <- t(t(x)*1e6/colSums(x))
  rpkm <- cpm * 1000 / len
  return(rpkm)
}

normalize <- function(raw_counts, gene_lengths, method='tpm', duplicates='combine',
                      verbose=TRUE, try_name_update=TRUE, out_dir_name='output'){
  initial_number_of_genes <- nrow(raw_counts)
  number_of_gene_lengths <- nrow(gene_lengths)
  
  gene_lengths$median_length <- as.integer(gene_lengths$median_length)
  
  # See how many of the expressed genes have known transcript lengths.
  expressed_gene_names <- raw_counts$gene_name
  know_transcript_length_names <- gene_lengths$gene_name
  
  number_known_lengths <- sum(expressed_gene_names %in% know_transcript_length_names)
  
  # Try updating gene names and see if that increases the number of identified 
  # transcript lengths.
  if (try_name_update){
    updated_expressed_gene_names <- unique(update_symbols(raw_counts$gene_name))
    updated_transcript_length_names <- unique(update_symbols(raw_counts$gene_name))
    
    updated_number_known <- sum(updated_expressed_gene_names %in% updated_transcript_length_names)
  }
  
  updated <- FALSE
  
  if (try_name_update & updated_number_known > number_known_lengths){
    updated <- TRUE
    
    raw_counts$gene_name <- update_symbols(raw_counts$gene_name)
    gene_lengths$gene_name <- update_symbols(raw_counts$gene_name)
    
    duplicated_new_names <- expressed_gene_names[duplicated(raw_counts$gene_name), ]
    duplicated_old_names <- raw_counts$gene_name[duplicated(raw_counts$gene_name), ]
    
    # Remove duplicated all but the first instance of duplicated gene names.
    # Maybe remove one with lower counts instead?
    raw_counts <- raw_counts[! duplicated(raw_counts$gene_name), ]
    gene_lengths <- gene_lengths[! duplicated(gene_lengths$gene_name), ]
    old_names <- expressed_gene_names[! duplicated(raw_counts$gene_name), ]
    
    num_dup_expressed <- initial_number_of_genes - ncol(raw_counts)
    num_dup_lengths <- number_of_gene_lengths - ncol(gene_lengths)
    
    update_str <- paste0(c(
      "Updating gene names was found to increase the number or transcripts
      for which there exists a known length.\n\n\t", s(number_known_lengths),
      '->', s(updated_number_known), '; +', 
      s(updated_number_known - number_known_lengths), '\n\n', 
      "After updating,", s(num_dup_expressed), 'genes were found to be 
      duplicates, and were removed accordingly. See duplicates file.'),
      collapse =' ')
    
    if (verbose){
      print(update_str)
    }
    
    write.csv2(data.frame(old_names=duplicated_old_names, 
                          new_names=duplicated_new_names),
               file=paste0(c(out_dir_name, '/duplicates.csv')))
  }
  
  # Remove genes that are missing gene lengths.
  missing <- ! raw_counts$gene_name %in% gene_lengths$gene_name
  missing_genes <- raw_counts$gene_name[missing]
  if (verbose){
    print(paste0('Omitting', as.character(sum(missing)), 
                 'genes, as their lengths were not found. See missing_genes.csv
                 for list of genes that were ommitted due to'), 
          collapse = ' ')
  }
  
  raw_counts <- raw_counts[!missing, ]
  
  # Save all genes for which lengths were not found to a CSV file. These genes
  # are pruned from the analysis.
  if (! updated){
    write.csv(data.frame(missing_genes=missing_genes), 
              file='output/missing_genes.csv')
  } else {
    write.csv(data.frame(new_names=missing_genes, 
                         old_names=old_names[missing]), 
              file='output/missing_genes.csv')
  }
  
  # Remove gene lengths that we will not need.
  missing <- ! gene_lengths$gene_name %in% raw_counts$gene_name
  missing_genes <- gene_lengths$gene_name[missing]
  gene_lengths <- gene_lengths[!missing, ]
  
  # Align gene lengths names to reads so that they are directly comparable.
  map <- match(gene_lengths$gene_name, raw_counts$gene_name)
  raw_counts <- raw_counts[map,]
  
}
  
  
  
  
  
  


# Check that they are properly aligned.
which(gene_lengths$gene_name != raw_counts$gene_name)

```

# Get TPMs
```{r}

tpms <- raw_counts
tpms[,2:ncol(tpms)] <- tpm3(raw_counts[,2:ncol(tpms)], gene_lengths[,2])
write.table(data.frame(tpms=tpms), file='output/tpms.tsv',
            sep='\t')

unique(as.integer(apply(tpms[,2:length(tpms)], 2, sum)))

```


