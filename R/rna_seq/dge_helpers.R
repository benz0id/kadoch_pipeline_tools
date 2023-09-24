library(ggplot2)


#' Convert Multiple Files Containing Read Counts into a Single Read Counts DF
#'
#' @param data_filepath The file containing the reads files.
#' 
#' @param pattern The pattern found in each of the reads files.
#' 
#' @param data_column The column form which to eaxtract the raw counts in each
#' of the reads files.
#' 
#' @param num_desc_col The number of columns that describe features of each gene.
#' 
#' @param gene_name_col The column that contains gene names.
#' 
#' @param gene_name_start The row at which gene names (and their reads) start.
#'
#' @return A dataframe containing the raw counts from each of the raw counts
#' files.
#' 
#' @export
#'
#' @examples
read_counts_files <- function(data_filepath, pattern, data_column, num_desc_col,
                              gene_name_col, gene_name_start){
  
  data_filepath <- paste(c(data_filepath, '/'), collapse='')
  
  # Load in data.
  samples_names <- list.files(data_filepath, pattern = pattern)
  
  # Construct storage dataframe for samples.
  sample_one <- readr::read_delim(paste0(data_filepath, samples_names[1]), 
                                  delim = "\t", escape_double = FALSE, 
                                  trim_ws = TRUE, col_names = FALSE, show_col_types = FALSE, progress=FALSE)
  gene_names <- sample_one[GENE_NAME_START:length(unlist(sample_one[,1])), gene_name_col]
  
  raw_counts <- data.frame(gene_name = gene_names)
  colnames(raw_counts) <- c('gene_name')
  
  # Parse counts into samples data frame.
  for (i in seq_along(samples_names)){
    sample_filename <- samples_names[i]
    
    # Extract counts from sample.
    sample_data <- readr::read_delim(paste0(data_filepath, sample_filename), 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE, col_names = FALSE, show_col_types = FALSE)
    sample_counts <- sample_data[GENE_NAME_START:length(unlist(sample_one[,1])), DATA_COLUMN]
    
    # Extract sample ID from sample.
    sample_ID <- unlist(strsplit(sample_filename, split = '-'))[1]
    
    raw_counts[,i + num_desc_col] <- sample_counts
    colnames(raw_counts)[i + num_desc_col] <- sample_ID
    
  }
  return(raw_counts)
}


#' Fetch Sample Treatment Conditions From a Sample Sheet Given the Sample IDs.
#'
#' @param sample_sheet_path Path to the sample sheet.
#' @param sample_id The id of the given sample.
#' @param names_start_row The row at which the first name can be found.
#' @param data_start The starting position of the treatment condition, relative
#' to the start of the sample name.
#' @param data_end The end position of the treatment condition, in units from
#' the end of the sample name.
#' @param out_delim How to separate the condition components in the output.
#'
#' @return A string representing the treatment that was found.
#' @export
#'
#' @examples
id_to_condition <- function(sample_sheet_path, sample_id, names_start_row, 
                            data_start, data_end, out_delim='&'){
  
  RNA_SampleSheet <- readr::read_csv(sample_sheet_path, show_col_types = FALSE,
                                     progress = FALSE)
  
  # Ensure SAMPLE_SHEET_NAMES_START is correct.
  if (!grepl("Sample_ID", RNA_SampleSheet[(names_start_row - 1), 1])){
    stop("SAMPLE_SHEET_NAMES_START is incorrectly configured.")
  }
  
  # Find coresponding complete sample ID
  sample_ids <- unlist(RNA_SampleSheet[names_start_row:length(unlist(RNA_SampleSheet[1])), 1])
  pos <- grep(sample_id, sample_ids)
  sample_name <- sample_ids[pos]
  sample_id_components <- unlist(strsplit(sample_name, split='_')) 
  condition_components <- sample_id_components[data_start:(length(sample_id_components) - data_end)]
  condition <- paste(condition_components, collapse=out_delim)
  
  return(condition)
}


#' Extract the Treatment Parameters from a Sample Name
#'
#' @param sample_name The full sample name.
#' @param expected_components The number of expected components, excess 
#' components will be considered to be treatment parameters.
#'
#' @return
#' @export
#'
#' @examples
extract_condition <- function(sample_name, expected_components){
  sample_id_components <- unlist(strsplit(sample_name, split='_')) 
  condition_components <- sample_id_components[(expected_components - 1):(length(sample_id_components) - 1)]
  condition <- paste(condition_components, collapse='&')
  return(condition)
}

#' Get a Bar Graph Displaying the Counts of a Specific Gene Across Conditions
#'
#' @param expression_data The expression data of interest, can be RPKM, CPM or 
#' even raw counts.
#' @param gene_name The name of the gene of interest.
#' @param sample_descs Slightly more detailed description of each sample, to be
#' displayed on the Y axis.
#' @param groups The group to which sample belongs.
#' @param title The title for the graph.
#' @param x The title for the x axis.
#' @param y The title for the y axis.
#' @param fill The title for the fill bar.
#' @param num_desc_cols The number of descriptive columns that prefix the 
#' <expression_data>.
#'
#' @return A horizontal barplot.
#' @export
#'
#' @examples
get_gene_counts_bar <- function(expression_data, gene_name, sample_descs,
                                groups, title, x, y, fill="Condition", 
                                num_desc_cols=1, colours=NULL){
  gene_row <- which(gene_name == expression_data$gene_name)
  
  if (length(gene_row) == 0){
    stop(paste0(c(gene_name, 'not found in expression data.'), collapse = ' '))
  } else if (length(gene_row) != 1){
    stop(paste0(c("Multiple matches found for", gene_name, "in expression data:",
                  expression_data$gene_name[gene_row]), collapse = ' '))
  }
  
  counts <- expression_data[gene_row, (num_desc_cols + 1): ncol(expression_data)]
  
  plt <- get_basic_bar(as.numeric(counts), sample_descs, groups, title, x, y, fill, colours)
  
  return(plt)
}

#' Generate A Basic Bar Plot
#'
#' @param values Some number of integer values.
#' @param value_labels Labels for each value.
#' @param groups The group to which each value belongs.
#' @param title The title of the barplot
#' @param x The X axis title.
#' @param y The Y axis title
#' @param fill The title for the fill bar.
#'
#' @return A horizontal barplot.
#' @export
#'
#' @examples
get_basic_bar <- function(values, value_labels, groups, title, x, y,
fill, colours=NULL){
  if (is.null(colours)){
    n <- length(unique(groups))
    colours <- colorspace::heat_hcl(n, 100, l=c(50, 90), power=1)
  }
  
  if (! is.factor(value_labels)){
    value_labels <- factor(value_labels, ordered = TRUE,
                           levels=rev(value_labels))
  }
  
  if (! is.factor(groups)){
    groups <- factor(groups, ordered = TRUE,
                           levels=unique(groups))
  }
  
  counts_df <- data.frame(label=value_labels, bars=values, group=groups)

  plot <- ggplot(counts_df, aes(x = label, y = bars, fill = group)) +
    scale_y_continuous(labels = scales::label_number(), limits = c(0, max(values))) +
    geom_col() +
    labs(title = title,
         x = x, y = y, fill = fill) +
    theme(plot.title = element_text(face = "bold")) +
    coord_flip() + 
    scale_fill_manual(values=colours) +
    theme_bw()
  
  
  return(plot)
}



#' Reorder Columns Belonging to Groups
#'
#' @param dataframe A dataframe containing n columns.
#' 
#' @param groups A vector of categories to which each column belongs,
#'  i.e. dataframe[, i] represents a sample belonging to groups[i].
#'    length(groups) == ncol(dataframe)
#'  
#' @param group_order A vector representing the order of the values in groups.
#'  where 
#'    length(group_order) <= groups
#'    all(group_order %in% unique(groups))
#'
#' @return A dataframe ordered by group_order
#' @export
#'
#' @examples
reorder_df_columns <- function(dataframe, groups, group_order, num_desc_columns=1){
  fac <- factor(groups, ordered=TRUE,levels = group_order)
  return(cbind(dataframe[1], dataframe[, (order(fac) + 1)]))
}


#' Reorder A Vector
#' 
#' @param groups A vector of categories (any type).
#'  
#' @param group_order A vector representing the order of the values in groups.
#'  where 
#'    length(group_order) <= groups
#'    all(group_order %in% unique(groups))
#'
#' @return A dataframe ordered by group_order
#' @export
#'
#' @examples
reorder_groups <- function(groups, group_order){
  fac <- factor(groups, levels = group_order, ordered=TRUE)
  return(groups[order(fac)])
}


#' Perform Basic Preliminary Analyses on Read Count Information
#'
#' @param raw_counts Raw count data.
#' @param groups The groups to which those counts belong.
#'
#' @return
#' @export
#'
#' @examples
reads_analysis <- function(raw_counts, groups, num_desc_columns=1, colours=NULL){
  
  # Extract count information.
  counts_per_samples <- apply(raw_counts[2:ncol(raw_counts)], MARGIN=2, FUN = sum)
  print(colours)
  plt <- get_basic_bar(values = counts_per_samples, names(counts_per_samples),
                              groups, title = "Total Number of Reads in Each Sample",
                              x = "Sample", y = "Total Reads", fill = "Treatment",
                              colours = colours)
  return(plt)
}

#' Normalize Counts Using the TMM Method.
#'
#' @param raw_counts Raw counts data, formatted as a dataframe.
#' @param groups The treatment groups, one for each sample in raw_counts.
#' @param show_factos Display a tibble of normalisation factors.
#' 
#' @return Normalised counts, in units of log2 cpm.
#'
#' @examples
tmm_normalize <- function(raw_counts, groups, num_desc_columns, show_factors=TRUE){
  # Convert to matrix & calculate normalization factors.
  filtered_data_matrix <- as.matrix(raw_counts[,(num_desc_columns + 1):ncol(raw_counts)])
  rownames(filtered_data_matrix) <- raw_counts$gene_name
  d <-  edgeR::DGEList(counts=filtered_data_matrix, group=groups)
  normalized <-  edgeR::calcNormFactors(d, method='TMM')
  
  if(show_factors){
    print(data.frame(Sample=colnames(filtered_data_matrix),
                     Norm_Factor=normalized$samples$norm.factors))
  }
  
  # Apply normalization factors and convert to log2 cpms.
  normalized <- limma::voom(normalized)
  norm <- cbind(raw_counts[1], normalized)
  return(norm)
}




#' Quantify Differential Expression From a Formatted Data Frame
#' 
#' Uses linear modelling provided by limma to compute the likelehood of 
#'
#' @param norm_counts A dataframe containing normalized read counts for the 
#' samples of interest (log2 cpm). Must have a column named "gene_name" 
#' containing the names of each gene in the data frame
#' 
#' @param compare_groups The groups to be compared. Should contain a name for
#' each sample in <norm_counts>. Considered to be a model of the data.
#' 
#' @param num_desc_columns The number of descriptive columns that lead the 
#' <norm_counts>.
#' 
#' @param sort_by_pval Whether to sort the output by significance.
#'
#' @return A list of genes ranked by adjusted significance of differential
#' expression.
#' 
#'
#' @examples
#'
#'
#'
get_de_lm <- function(norm_counts, compare_groups, num_desc_columns, 
                      sort_by_pval=FALSE){
  # Use linear models from limma package to identify differentially expressed genes.
  
  # norm_counts <- norm_cpms
  # sample_groups <- samples
  
  # Define our model.
  groups_to_compare <- unique(compare_groups)
  
  if (length(groups_to_compare) != 2){
    stop("Two unique groups are required.")
  }
  
  fg <- compare_groups == groups_to_compare[1]
  
  model_design <- as.matrix(data.frame(intercept=1, sg=!fg))
  
  # Format expression data before fitting.
  expressionMatrix <- as.matrix(norm_counts[, (num_desc_columns + 1):ncol(norm_counts)])
  rownames(expressionMatrix) <- norm_counts$gene_name
  colnames(expressionMatrix) <- colnames(norm_counts)[(num_desc_columns + 1):length(norm_counts)]
  minimalSet <- Biobase::ExpressionSet(assayData=expressionMatrix)
  
  # Generate linear model.
  fit <- limma::lmFit(minimalSet, model_design)
  
  # Apply eBayes to compute probability of observing given data given the null.
  fit2 <- limma::eBayes(fit,trend=TRUE)
  
  # Convert to table.
  topfit <- limma::topTable(fit2, 
                            coef=ncol(model_design),
                            adjust.method = "BH",
                            number = nrow(expressionMatrix))
  
  
  # Format results; sort by pvalue.
  output_hits <- merge(norm_counts[,1:num_desc_columns],
                       topfit,
                       by.y=0,by.x=1,
                       all.y=TRUE)
  
  if (sort_by_pval){
    output_hits <- output_hits[order(output_hits$adj.P.Val, decreasing = FALSE),]
  }
  
  colnames(output_hits)[1] <- 'gene_name'
  
  return(output_hits)
}


#' Quantify Differential Expression Between Two Conditions
#' 
#' Returns the significance of differential expression for each gene between
#' the two given conditions. Fold-change values are relative to <reference_cond>,
#' that is, a gene that is expressed twice as much in compare_cond will have
#' a reported logFC value of 1.
#'
#' @param expression_data A dataframe containing normalized read counts for the 
#' samples of interest (log2 cpm). Must have a column named "gene_name" 
#' containing the names of each gene in the data frame.
#' 
#' @param reference_cond The string name of the first condition to be compared. 
#' Groups must contain at least one element equal to <reference_cond>.
#' 
#' @param compare_cond The string name of the first condition to be compared. 
#' Groups must contain at least one element equal to <reference_cond>.
#' 
#' 
#' @param groups The a vector of treatments/conditions, containing on entry for 
#' each sample in the dataframe.
#' 
#' @param col_offset The number of descriptive columns that prefix 
#' 
#' @param save_output 
#'
#' @return
#' @export
#'
#' @examples
get_de_between_conditions <- function(expression_data, reference_cond, 
                                      compare_cond, groups,
                                      col_offset=1, save_output=TRUE){
  # Where expression_data is a dataframe of genes CPMs and c1_name and c2_name
  # are names in group. col_offset is the number of leading description columns 
  # in <expression_data>.
  
  # Extract and format expression dataset and groupings.
  
  formatted_data <- cbind(expression_data[1:col_offset], 
                          get_cols(expression_data, groups,reference_cond, col_offset),
                          get_cols(expression_data, groups,compare_cond, col_offset))
  
  colnames(formatted_data) <- c('gene_names', letters[1:(ncol(formatted_data) - 1)])
  
  formatted_groups <- rep(reference_cond, sum(groups == reference_cond))
  formatted_groups <- append(formatted_groups, rep(compare_cond, 
                                                   sum(groups == compare_cond)))
  
  
  output_hits <- get_de_lm(formatted_data, formatted_groups, col_offset)
  
  if (save_output){
    out_name <- paste0(c('DE_list_', reference_cond, '_', compare_cond, '.csv'), 
                       collapse='')
    write.csv(output_hits, paste0('output/', out_name))
  }
  
  return(output_hits)
}


#' Get a Boolean Vector of DEGs.
#'
#' @param expression_data A dataframe containing normalized read counts for the 
#' samples of interest (log2 cpm). Must have a column named "gene_name" 
#' containing the names of each gene in the matrix.
#' 
#' @param reference_cond The string name of the first condition to be compared.
#' 
#' @param compare_cond The string name of the second condition to be compared.
#' 
#' @param groups The a vector of treatments/conditions, containing on entry for 
#' each sample in the dataframe.
#' 
#' @param adjpval_threshold The adjusted p-value threshold for inclusion into 
#' the list of DE genes
#' 
#' @param logfc_threshold The log fold change threshold for inclusion into the
#' list of DE genes. Must be a positive number. Genes with sufficient FC in 
#' either direction will be returned.
#' 
#' @param thr_set One of "increase", "decrease", or "abs". Returns genes with
#' a positive fold change, negative fold change, or absolute with <thr_set> 
#' magnitude.
#'
#' @return A binary named vector, where each entry is a boolean representing
#' whether that gene is DE.
#' 
#' @export
#'
#' @examples
get_de_binary <- function(expression_data, reference_cond, compare_cond, groups, 
                          adjpval_threshold, logfc_threshold, thr_set='abs'){
  # Get a binary named vector, where each entry stores whether each gene is DE,
  # described above.
  de <- get_de_between_conditions(expression_data, reference_cond, compare_cond,
                                  groups)
  de_gene_names <- get_de_gene_names(de, adjpval_threshold, logfc_threshold,
                                     thr_set = thr_set)
  
  is_de_expressed <- expression_data$gene_name %in% de_gene_names
  names(is_de_expressed) <- unlist(expression_data$gene_name)
  return(is_de_expressed)
  
}


#' Return the Genes that Lie at the Exclusive Intersection of number of DE 
#'
#' @param expression_data Filtered and normalized expression data.
#' @param controls A list of n control conditions - corresponding to each 
#' treatment.
#' @param treatments A list of n treatment conditions.
#'
#' @return
#' @export
#'
#' @examples
get_venn_df <- function(expression_data, controls, treatments, groups, 
                        adjpval_threshold, logfc_threshold, thr_set='abs',
                        include_control_name=FALSE, out_file=NA){
  # Iterate through all subsets of the diagram.
  
  # Binary counter with least sig dig first used to iterate through all subsets 
  # of the diagram.
  counter <- logical(length(controls))
  
  # Used to iterate the counter.
  update_counter <- function(counter){
    j <- 1
    while (counter[j] == TRUE & j <= length(counter)){
      counter[j] <- FALSE
      j = j + 1
    }
    
    # Return if counter is complete.
    if (j > length(counter)){return(counter)} 
    else {counter[j] <- TRUE; return(counter)}
  }
  
  # Get the name of the current set, as denoted by the current state of
  # the counter.
  get_cur_intersection_name <- function(){
    get_set_name <- function(i){
      if (include_control_name){
        set_name <- paste0(c(controls[i], '_vs_', treatments[i]), collapse='')
      } else {
        set_name <- treatments[i]
      }
      return(set_name)
    }
    
    condition_inclusions <- c()
    for (i in seq_along(counter)){
      included <- counter[i]
      if (included){
        condition_inclusions <- c(condition_inclusions, get_set_name(i))
      } else {
        s <- paste0(c('not', get_set_name(i)), collapse=' ')
        condition_inclusions <- c(condition_inclusions, s)
      }
    }
    
    return(paste0(condition_inclusions, collapse = ' and '))
  }
  
  # Construct matrix of DE that tells us whether each gene is DE in each
  # condition.
  de_binary_inclusions <- matrix(ncol=length(treatments),
                                 nrow=nrow(expression_data))
  rownames(de_binary_inclusions) <- expression_data$gene_name
  for (i in seq_along(controls)){
    control <- controls[i]
    treat <- treatments[i]
    
    deb <- get_de_binary(expression_data, control,  treat, groups, 
                         adjpval_threshold, logfc_threshold, thr_set)
    
    if (! all(names(deb) == rownames(de_binary_inclusions))){
      stop("Misalignment between rownames when fetching DE genes.")
    }
    
    de_binary_inclusions[, i] <- deb
  }
  
  # Get all of the genes that are DE in all conditions currently denoted by 
  # counter.
  get_cur_intersection <- function(){
    incl <- as.logical(rep(TRUE, nrow(expression_data)))
    for (i in seq_along(treatments)){
      if (counter[i]){
        incl <- incl & de_binary_inclusions[, i]
      } else {
        incl <- incl & !de_binary_inclusions[, i]
      }
    }
    return(incl)
  }
  
  subset_inclusion_data <- data.frame(gene_name=expression_data$gene_name)
  j <- 1
  # While the counter has not completed, iterate through every permutation.
  while (! all(counter)){
    counter <- update_counter(counter)
    j <- j + 1
    
    subset_inclusion_data[, j] <- get_cur_intersection()
    names(subset_inclusion_data)[j] <- get_cur_intersection_name()
  }
  
  
  # Make sure that union of all sets is what we expect it to be.
  sub_matrix <- as.matrix(subset_inclusion_data[, 2:ncol(subset_inclusion_data)])
  sub_sums <- apply(sub_matrix, 1, sum)
  if (any(sub_sums) > 1){
    stop('Intersections calculated improperly.')
  }
  
  sub_ag <- apply(sub_matrix, 1, any)
  binary_ag <- apply(de_binary_inclusions, 1, any)
  if (any(sub_ag != binary_ag)){
    stop("Intersections calculated improperly.")
  }
  
  
  # Convert binary inclusion sets into gene lists.
  venn_genes <- data.frame(gene_name=subset_inclusion_data$gene_name)
  lens <- c()
  
  for (col in seq(2, ncol(subset_inclusion_data))){
    genes <- venn_genes$gene_name
    sel <- unlist(subset_inclusion_data[col])
    len <- length(genes[sel])
    lens <- c(lens, len)
    venn_genes[col] <- c(genes[sel], rep('', length(genes) - length(genes[sel])))
    
    
  }
  
  colnames(venn_genes) <- colnames(subset_inclusion_data)
  
  venn_genes <- venn_genes[-(max(lens):nrow(venn_genes)), -1]
  
  
  if (!is.na(out_file)){
    filename <- paste0(c(treatments, thr_set), collapse = '_')
    write.csv(venn_genes, paste0(c(out_file, '/', 'subset_info_', 
                                              filename, '.csv'),
                                            collapse = ''))
  }
  
  return(venn_genes)
}


#' Get the Gene Symbols of All DE Genes.
#'
#' @param de_gene_list A data frame containing all differentially expressed 
#' genes.
#' 
#' @param adjpval_threshold The adjusted p-value threshold for inclusion into 
#' the list of DE genes
#' 
#' @param logfc_threshold The log fold change threshold for inclusion into the
#' list of DE genes. Must be a positive number.
#' 
#' @param thr_set One of "increase", "decrease", or "abs". Returns genes with
#' a positive fold change, negative fold change, or absolute with <thr_set> 
#' magnitude.
#' 
#' @param bool 
#'
#' @return
#' @export
#'
#' @examples
get_de_gene_names <- function(de_gene_list, adjpval_threshold, logfc_threshold, 
                              thr_set='abs', bool=FALSE){
  # Return a vector containing all gene names with expression profiles matching
  # the given descriptions.
  # thr_set is one of "increase", "decrease", or "abs".
  
  if (logfc_threshold < 0){
    stop("change threshold must be absoute. specify increase or decrease 
         using the <thr_set> argument.")
  }
  
  sig <- de_gene_list$adj.P.Val <= adjpval_threshold
  
  gt <- de_gene_list$logFC >= logfc_threshold
  lt <- de_gene_list$logFC <= -logfc_threshold
  
  if (thr_set == 'abs'){
    valid <- sig & (gt | lt) & !(gt & lt)
  } else if (thr_set == 'increase') {
    valid <- sig & gt
  } else if (thr_set == 'decrease') {
    valid <- sig & lt
  } else {
    stop(paste0(c(thr_set, ' is not a valid threshold setting'), collapse = ''))
  }
  
  if (bool){
    names(valid) <- de_gene_list$x
    return(valid)
  }
  
  rtrn <- unlist(de_gene_list$gene_name)[valid]
  names(rtrn) <- c()
  
  return(rtrn)
  
}



DEFAULT_VOLCANO_TITLE <- "Identifying Key DE Genes Through Significance and Magnitude of DE"

#' Create a Volcano Plot
#'
#' @param de_gene_list A data frame containing all differentially expressed 
#' genes.
#' 
#' @param adjpval_threshold The adjusted p-value threshold for inclusion into 
#' the list of DE genes
#' 
#' @param logfc_threshold The log fold change threshold for inclusion into the
#' list of DE genes. Must be a positive number.
#' 
#' @param reference_cond The group/condition DE expression is relative to.
#' 
#' @param num_names The number of gene names to display on the the volcano plot.
#' 
#' @param title The title to display.
#' 
#'
#' @return ggplot of the volcano plot.
#' @export
#'
#' @examples
volcano <- function(de_gene_list, logfc_threshold, sig_threshold, reference_cond,
                    num_names=5, title=DEFAULT_VOLCANO_TITLE, 
                    genes_to_show=NULL, save_plot=TRUE){
  
  # Copy to avoid aliasing.
  de_gene_data <- data.frame(de_gene_list)
  
  # Convert to absolute threshold.
  logfc_threshold <- abs(logfc_threshold)
  
  # Which genes are significantly differentially expressed and have high abs fold change?
  de_gene_data$diffexpressed <- "NO"
  up <- get_de_gene_names(de_gene_data, sig_threshold, logfc_threshold, thr_set='increase', bool=TRUE)
  down <- get_de_gene_names(de_gene_data, sig_threshold, logfc_threshold, thr_set = 'decrease', bool=TRUE)
  de_gene_data$diffexpressed[up] <- "UP"
  de_gene_data$diffexpressed[down] <- "DOWN"
  
  #  Include <num_names> most "outlying" genes names in the plot.
  de_gene_data$delabel <- NA
  de_gene_data$fac <- abs(as.numeric(de_gene_data$logFC)) * - log2(as.numeric(de_gene_data$adj.P.Val))
  outlier_rank <- order(de_gene_data$fac, decreasing = TRUE)
  most_interesting <- outlier_rank[1:num_names]
  de_gene_data$delabel[most_interesting] <- de_gene_data$gene_name[most_interesting]
  
  # Include requested genes in the plot.
  if (! is.null(genes_to_show)){
    if (! all(genes_to_show %in% de_gene_data$gene_name)){
      stop(paste0(c(genes_to_show[! (genes_to_show %in% de_gene_data$gene_name)], 
                    "not found in the list of DE genes."), 
                  collapse =' '))
    }
    to_highlight <- de_gene_data$gene_name %in% genes_to_show
    de_gene_data$delabel[to_highlight] <- de_gene_data$gene_name[to_highlight]
  }
  
  # Generate a caption. Ignore NA genes.
  caption <- paste0(c('Fold change relative to ', reference_cond, 
                      ' | logFC Threshold = +/-', as.character(logfc_threshold), 
                      ' | Significance Threshold = ', as.character(sig_threshold),
                      ' | N Up = ', as.character(sum(up & ! is.na(up))),
                      ' | N Down = ', as.character(sum(down& ! is.na(down)))), 
                    collapse='')
  
  # Plot the data.
  plt <- ggplot(data=de_gene_data, aes(x=logFC, y=-log2(adj.P.Val), col=diffexpressed, label = delabel)) + 
    geom_point() + 
    geom_vline(xintercept=c(-logfc_threshold, logfc_threshold), col="red") +
    geom_hline(yintercept=-log2(sig_threshold), col="red") +
    scale_color_manual(values=c('blue','grey', 'red')) +
    geom_label(na.rm=TRUE) +
    labs(title=title,
         caption=caption,) + 
    theme_bw() +
    theme(legend.position = 'none')
  
  if (save_plot){
    ggsave(paste0(c('output/', title, '_volcano_plot.svg'), collapse = ''),
           plot=plt, width=5, height=5)
  }
  
  return(plt)
}



#' Generate a Simple Title for Two conditions
#'
#' @param reference_cond String representation of the first conditon.
#' @param compare_cond String representation of the second conditon.
#'
#' @return
#'
#' @examples
get_title <- function(reference_cond, compare_cond){
  return(paste0(c(reference_cond, ' vs. ',  compare_cond), collapse=''))
}


#' Generate a Simple Title for Two conditions
#'
#' @param cond_names Vector of condition names.
#'
#' @return
#'
#' @examples
get_title_vec <- function(cond_names){
  if (length(cond_names) == 2){
    s <- ' and '
  } else{
    s <- ', and '
  }
  frst <- paste0(cond_names[seq(1, length(cond_names) - 1)], collapse=', ')
  return(paste0(c(frst, cond_names[length(cond_names)]), collapse=s))
}

#' Visualize Density of Read Lengths
#'
#' @param norm_counts A dataframe containing normalized read counts for the 
#' samples of interest (log2 cpm). Must have a column named "gene_name" 
#' containing the names of each gene in the data frame.

#' @param legend Boolean, whether to display a legend.

#' @param title Title of the plot to be displayed.
#' 
#' @param xlab The lab to place on the x axis
#'
#' @return nothing
#' @export
#'
#' @examples
get_count_density <- function(data2plot, legend=TRUE, title = '', 
                              xlab="log2-CPM"){
  # Function to help visualize the distribution of RNA lenghts in our dataset.
  
  # Display histogram of count densities.
  counts_density <- apply(data2plot, 2, density)
  
  # Calculate the limits across all the samples
  xlim <- 0; ylim <- 0
  for (i in 1:length(counts_density)) {
    xlim <- range(c(xlim, counts_density[[i]]$x));
    ylim <- range(c(ylim, counts_density[[i]]$y))
  }
  cols <- rainbow(length(counts_density))
  ltys <- rep(1, length(counts_density))
  
  # Plot the first density plot to initialize the plot
  p <- plot(counts_density[[1]], xlim=xlim, ylim=ylim, type="n",
            ylab="Density", xlab=xlab, main = title, cex.lab = 0.85)
  
  
  #plot each line
  for (i in 1:length(counts_density)){
    lines(counts_density[[i]], col=cols[i], lty=ltys[i])
  }
  #create legend
  if (legend){
    legend("topright", colnames(data2plot),
           col=cols, lty=ltys, cex=0.75,
           border ="blue", text.col = "green4",
           merge = TRUE, bg = "gray90")
  }
  
  return(p)
  
}

#' Compare Differential Expression Between all Conditions.
#'
#' @param norm_counts A dataframe containing normalized read counts for the 
#' samples of interest (log2 cpm). Must have a column named "gene_name" 
#' containing the names of each gene in the data frame.
#' 
#' @param groups The a vector of treatments/conditions, containing on entry for 
#' each sample in the dataframe.
#' 
#' @param adjpval_threshold The adjusted p-value threshold for inclusion into 
#' the list of DE genes
#' 
#' @param logfc_threshold The log fold change threshold for inclusion into the
#' list of DE genes. Must be a positive number.
#' 
#' @param num_names The number of gene names to display on the the volcano plot.
#'
#' @return a list containing all plots created.
#' @export
#'
#' @examples
compare_all <- function(norm_counts, groups, logfc_threshold, sig_threshold, num_names=4){
  unique_groups <- unique(groups)
  num_groups <- length(unique_groups)
  
  g1s <- c()
  g2s <- c()
  de_data <- list()
  de_gene_lists <- list()
  de_plots <- list()
  
  for (i in seq(1, length(unique_groups) - 1)){
    for (j in seq((i + 1), length(unique_groups))){
      
      reference_condition <- unique_groups[i]
      cond2 <- unique_groups[j]
      
      title <- get_title(reference_condition, cond2)
      
      g1s <- append(g1s, reference_condition)
      g2s <- append(g2s, cond2)
      
      de_datum <- get_de_between_conditions(norm_counts, reference_condition, cond2, groups)
      de_data <- rlist::list.append(de_data, de_datum)
      
      png(file=paste0(c('output/', title, '_volcano_plot.png'), collapse = ''),
          width = 1500, height = 700 )
      de_plot <- volcano(de_datum, logfc_threshold, sig_threshold, 
                         reference_cond=reference_condition, title=title, 
                         num_names=num_names)
      print(de_plot)
      dev.off()
      de_plots[[`title`]] <- de_plot
    }
  }
  
  return(de_plots)
}


#' Create Venn Diagrams Comparing DE Between Conditions
#'
#' @param norm_counts A dataframe containing normalized read counts for the 
#' samples of interest (log2 cpm). Must have a column named "gene_name" 
#' containing the names of each gene in the data frame.
#' 
#' @param control_group The string name of the control group in <groups>. Accepts
#' multiple control groups, in which case each control is assigned to the the
#' treatment at its respective index.
#' 
#' @param treatment_groups A vector containing the string names of each of the 
#' groups to be compared.
#' 
#' @param adjpval_threshold The adjusted p-value threshold for inclusion into 
#' the list of DE genes
#' 
#' @param logfc_threshold The log fold change threshold for inclusion into the
#' list of DE genes. Must be a positive number.
#' 
#' @param groups The a vector of treatments/conditions, containing on entry for 
#' each sample in the dataframe.
#' 
#' @param thr_set One of "increase", "decrease", or "abs". Returns genes with
#' a positive fold change, negative fold change, or absolute with <thr_set>
#' magnitude.
#'
#' @return the venn diagram object.
#' @export
#'
#' @examples
get_de_venn <- function(norm_cpms, control_group, treatment_groups, pval_threshold,
                        logfc_threshold, groups, thr_set, height=1000, width=1000,
                        margin=0.05, res=300, show_controls=FALSE, do_main=TRUE){
  # Generate a Venn diagram that compares the genes in common between each of
  # the treatment groups.
  
  control_vs_treatment <- list()
  
  cntrl <- control_group
  for (i in seq(1, length(treatment_groups))){
    if(length(control_group) != 1){
      cntrl <- control_group[i]
    }
    
    cond <- treatment_groups[i]
    title <- get_title(cntrl, cond)
    de_datum <- get_de_between_conditions(norm_cpms, cntrl, cond, groups)
    diff_reg <- get_de_gene_names(de_datum, adjpval_threshold, logfc_threshold, 
                                  thr_set=thr_set)
    control_vs_treatment[[`title`]] <- diff_reg
  }
  
  treatment_names_str <- paste0(treatment_groups, collapse = '_')
  filename <- paste0(c('output/', treatment_names_str, '_', 
                       thr_set,'_venn_diagram.png'), collapse = '')
  area_vec <- unlist(lapply(control_vs_treatment, length))
  fill <- c('lightpink', 'lightblue', 'lightgreen')[seq(1, length(control_vs_treatment))]
  
  if (thr_set == 'abs'){
    prefix <- 'Genes Differentially Expressed in'
  } else if (thr_set == 'increase') {
    prefix <- ' Genes Up Regulated in'
  } else if (thr_set == 'decrease') {
    prefix <- 'Genes Down Regulated in'
  } else {
    stop(paste0(c(thr_set, ' is not a valid threshold setting'), collapse = ''))
  }
  
  if (length(control_group) == 1){
    post <- paste0(c('Relative to', control_group), collapse=' ')
  } else {
    post <- ''
  }
  
  if(do_main){
    main <- paste0(c(prefix, get_title_vec(treatment_groups), post), collapse=' ')
  } else {
    main <- NULL
  }
  
  
  invisible(VennDiagram::venn.diagram(
    imagetype="tiff" ,
    height = height , 
    width = width, 
    resolution = res,
    compression = "lzw",
    x = control_vs_treatment,
    category.names = names(control_vs_treatment),
    filename = filename,
    output = TRUE,
    fill = fill,
    disable.logging = TRUE,
    ext.text = FALSE,
    margin=margin,
    main = main,
    disable_loggin=TRUE
  ))
  return(NULL)
}


#' Display all Images in a Directory with names Matching a Pattern
#'
#' @param directoy The target directory.
#' @param pattern All images with names matching <pattern> will be displayed.
#'
#' @return NULL
#' @export
#'
#' @examples
include_all_matching_images <- function(directoy, pattern){
  matching_images <- list.files(directoy, pattern = pattern)
  for(image_name in matching_images){
    image_path <- paste(c(directoy, '/', image_name), collapse='')
    knitr::include_graphics(image_path)
  }
}


#' Extract Columns From a Dataframe By Group
#'
#' @param targ_df The dataframe to extract columns from.
#' @param groups Group to which each column in the dataframe belongs.
#' @param group_name The name of the group of interest.
#' @param col_offset The number of columns leading the dataframe, containing 
#' gene identifiers.
#'
#' @return Dataframe composed of columns that belong to <group_name>.
#' @export
#'
#' @examples
get_cols <- function(targ_df, groups, group_name, col_offset){
  if (! any(groups == group_name)){
    stop(paste0(c("No groups matching \"", group_name, 
                  "\" found in the given groups."), collapse=''))
  }
  
  treatment_cols <- which(groups == group_name) + col_offset
  return(targ_df[treatment_cols])
}


# Convert a sample name to a sample ID.
name_to_sample_id <- function(sample_name){
  return(groups_df[1][groups_df[2] == sample_name])
}

# Convert a sample's name into it's treatment group.
name_to_group <- function(sample_name){
  components <- unlist(strsplit(sample_name, '_'))
  rep_info_loc <- grepl('Rep', components, ignore.case = TRUE)
  return(paste0(components[!rep_info_loc], collapse='_'))
}


#' Save an Image to The Output Folder
#'
#' @param image The plot, or any other image that can be displayed.
#' @param width Width of the image to be saved, in pixels.
#' @param height Height of the image to be saved, in pixels.
#' @param name Name of the image.
#' @param display Whether to display the image after saving.
#'
#' @return
#' @export
#'
#' @examples
save_img <- function(image, width, height, name, display=TRUE){
  png(file=paste0(c('output/', name), collapse = ''),
      width = width, height = height )
  print(image)
  dev.off()
  if(display){print(image)}
}
