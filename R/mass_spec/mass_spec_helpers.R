DEFAULT_FILE_START <- 4
DATA_START <- 2

neg_ind <- function(index, len){
  if (index < 0){
    return(len - index)
  } else {
    return(index)
  }
}

#' Get Marks from a Data Frame or CSV
#' 
#' @param data_frame
#' @param start The start position of the mark after separating by sep.
#' @param stop The stop position of the mark after separating by sep.
#' @return A character vector.
get_ms_marks <- function(data_frame, start, stop, sep='_'){
  start <- neg_ind(start)
  stop <- neg_ind(stop)
  
  if (is.character(data_frame)){
    data_start <- 3
    raw_data <- readr::read_csv(csv_path)
    columns <- colnames(raw_data)
    
    names <- c()
    
    i <- data_start
    while (grepl('Unique', columns[i]) | 
           grepl('Total', columns[i])){
      names <- c(names, columns[i])
      i <- i + 2
    }
  } else if (is.data.frame(data_frame)){
    names <- colnames(data_frame)[DATA_START:ncol(data_frame)]
  }
  
  marks <- c()
  for (name in names){
    name <- unlist(strsplit(name, ' '))[2]
    mark <- paste(unlist(strsplit(name, sep))[start:stop], collapse=sep)
    marks <- c(marks, mark)
  }
  return(marks)
}


#' Parse a Mass Spec DF
#' 
#' @param csv_path Path to the csv formatted matrix.
#' @param sel_unique Whether to return unique peptides.
#' @return The formatted data frame, with gene names in the first column.
#' 
parse_mass_spec_matrix <- function(csv_path, sel='total', sp_filter=TRUE,
                                   add_refs_for_missing=T, gene_name_col=2){
  data_start <- 5
  raw_data <- readr::read_csv(csv_path)
  columns <- colnames(raw_data)
  
  out_df <- data.frame(gene_name=raw_data[gene_name_col])
  
  # If no gene name is available, use reference.
  if (add_refs_for_missing){
    nas <- is.na(out_df$gene_name)
    out_df$gene_name[nas] <- raw_data$reference[nas]
  }
  
  # Find requested columns.
  data_end <- data_start
  if (sel != 'intensity'){
    data_start <- 5
    data_end <- data_start
    while (grepl('Unique', columns[data_end]) | 
           grepl('Total', columns[data_end])){
      data_end <- data_end +  1
    }
    data_end <- data_end - 1
  } else if (sel == 'intensity') {
    data_start <- 5
    while(! grepl('Sum Intensity', columns[data_start])){
      data_start <- data_start + 1
    }
    data_end <- data_start
    while(grepl('Sum Intensity', columns[data_end])){
      data_end <- data_end + 1
    }
  } else {
    stop('Unrecognised data type')
  }
  
  for(i in seq(data_start, data_end)){
    is_unique <- grepl('Unique', columns[i])
    
    if (sel == 'total' & grepl('Total', columns[i])){
        incl <- T
    } else if (sel == 'intensity' & grepl('Sum Intensity', columns[i])){
        incl <- T
    } else if (sel == 'unique' & grepl('Unique', columns[i])){
        incl <- T
    } else {
        incl <- F
    }
    
    print(paste(columns[i], is_unique))
    
    if (incl){
      out_df <- cbind(out_df, raw_data[,i])
      name <- unlist(strsplit(colnames(raw_data)[i], ' '))[2]
    }
  }
  
  if (sel == 'intensity'){
    for (i in seq(2, ncol(out_df))){
      for (j in seq(nrow(out_df))){
        if (out_df[j, i] == 'NF'){
          out_df[j, i] <- 0
        }
      }
    }
    for(col in seq(2, ncol(out_df))){
      out_df[col] <- as.numeric(out_df[[col]])
    }
  }
  
  if (sp_filter){
    keep <- grepl('sp\\|', raw_data$reference)
    print(paste(sum(keep), sum(!keep)))
    out_df <- out_df[keep,]
  }
  
  colnames(out_df)[1] <- 'gene_name'
  return(out_df)
}


#' Normalize by IgG in an MS Dataframe.
#' 
#' @param peptide_counts A data frame with genes in the first column and counts for 
#' each condition in each other or matrix with gene names as row names.
#' Must have an IgG for each treatment.
#' 
#' @param method How to normalize to IgG. Either subtraction or division.
#' 
#' @return A data frame of the same format as the input dataframe, without the
#'  input columns.
#'  
igg_normalize <- function(peptide_counts, marks, conditions, method='sub',
                          norm_mark='igg'){
  if (is.data.frame(peptide_counts)){
    counts_matrix <- as.matrix(peptide_counts[DATA_START:ncol(peptide_counts)])
    rownames(counts_matrix) <- peptide_counts$gene_name
  }
  else{
    counts_matrix <- peptide_counts
  }
  
  is_igg <- grepl(norm_mark, marks, ignore.case = T)
  print(is_igg)
  
  new_marks <- marks[! is_igg]
  new_conditions <- conditions[! is_igg]
  normalized_matrix <- counts_matrix[,! is_igg]
  
  for (i in seq_along(new_conditions)){
    for (j in seq_along(conditions)){
      cond_match <- new_conditions[i] == conditions[j]
      if (cond_match & is_igg[j]){
        if (method == 'div'){
          print(paste(new_marks[i], new_conditions[i], ' / ', marks[j], conditions[j]))
          normalized_matrix[,i] <- normalized_matrix[,i] / counts_matrix[,j]
        } else if (method == 'sub') {
          print(paste(new_marks[i], new_conditions[i], ' - ', marks[j], conditions[j]))
          normalized_matrix[,i] <- normalized_matrix[,i] - counts_matrix[,j]
        } else {
          stop('Unrecognised normalization method.')
        }
      }
    }
  }
  
  return(normalized_matrix)
}

#' Normalize a MS dataframe to the amount of target protein pulled down.
#' 
#' @param peptide_counts A data frame with genes in the first column and counts for 
#' each condition in each other or matrix with gene names as row names.
#' 
#' @return A data frame of the same format as the input dataframe, without the
#' input columns.
#'  
norm_to_target <- function(peptide_counts, marks, conditions, return_factors=F){
  if (is.data.frame(peptide_counts)){
    counts_matrix <- as.matrix(peptide_counts[DATA_START:ncol(peptide_counts)])
    rownames(counts_matrix) <- peptide_counts$gene_name
  }
  else{
    counts_matrix <- peptide_counts
  }
  
  norm_matrix <- matrix(data=0, 
                        nrow=nrow(counts_matrix), 
                        ncol=ncol(counts_matrix))
  
  norm_factors <- c()
  
  for (col in seq(ncol(counts_matrix))){
    mark <- marks[col]
    mark_row <- which(rownames(counts_matrix) == mark)
    if (length(mark_row) < 1){
      stop(paste(mark, 'not found in rows'))
    } else if (length(mark_row) > 1){
      stop(paste('Multiple rows found matching', mark))
    }
    
    norm_val <- counts_matrix[mark_row, col]
    norm_factors <- c(norm_factors, norm_val)
    norm_matrix[,col] <- counts_matrix[,col] / norm_val
  }
  
  
  if (return_factors){
    return(norm_factors)
  }
  print(norm_factors)
  
  colnames(norm_matrix) <- colnames(counts_matrix)
  rownames(norm_matrix) <- rownames(counts_matrix)
  
  return(norm_matrix)
}


#' Get the Rows Matching the Given list of Genes
#' 
#' @param peptide_counts A data frame with genes in the first column and counts for 
#' each condition in each other.
#' @param genes A vector of gene names/
#' @param as_matrix Whether to convert to a matrix before returning.
get_genes_data_frame <- function(peptide_counts, genes, as_matrix=F){
  
  if (is.matrix(peptide_counts)){
    data_frame <- cbind(rownames(peptide_counts), as.data.frame(peptide_counts))
    colnames(data_frame)[1] <- 'gene_name'
    rownames(data_frame) <- c()
  } else {
    data_frame <- peptide_counts
  }
  
  gene_inds <- c()
  for (gene in genes){
    ind <- which(data_frame$gene_name == gene)
    if (length(ind) != 1){
      print(paste('Unexpected number of genes found for', gene, length(ind), 
                  '!=', 1, '.', paste(data_frame$gene_name[ind]), '==', gene))
    }
    gene_inds <- c(gene_inds, ind)
  }
  new_data_frame <- data_frame[gene_inds,]
  print(new_data_frame$gene_name)
  if (as_matrix){
    new_matrix <- as.matrix(new_data_frame[,DATA_START:ncol(new_data_frame)])
    print(new_data_frame$gene_names)
    rownames(new_matrix) <- new_data_frame$gene_name
    return(new_matrix)
  }
  return(new_data_frame)
}

