#' Generate a .rnk File Using The Given DE List
#'
#' @param de_list 
#' @param ref_condition 
#' @param treat_codition 
#' @param out_dir 
#' @param return_dataframe If the specified rank file exists, return it, else
#' create a new one and return it without saving.
#'
#' @return A path to the rank file.
#' @export
#'
#' @examples
get_rank_list <- function(de_list, ref_condition, treat_condition, out_dir='output'){
  
  rankfile_path <- paste0(c(out_dir, '/', ref_condition, '-', treat_condition, '.rnk'), collapse='')
  
  rnk <- data.frame(gene_name=de_list$gene_name)
  
  dir_pvals <- numeric()
  for (i in seq_along(de_list$adj.P.Val)){
    logfc <- de_list$logFC[i]
    adjpval <- de_list$adj.P.Val[i]
    dir_pvals[i] <- -log10(adjpval) * sign(logfc)
  }
  
  de_list$rank <- dir_pvals
  
  rnk$dir.adj.P.Val <- dir_pvals
  
  rnk <- rnk[order(rnk$dir.adj.P.Val, decreasing=FALSE),]
  
  write.table(rnk, rankfile_path, sep='\t', row.names=FALSE, col.names = FALSE)
  
  return(rankfile_path)
}

#' Find a .rnk File Using The Given Attributes
#'
#' @param de_list 
#' @param ref_condition 
#' @param treat_codition 
#' @param out_dir 
#' @param return_dataframe If the specified rank file exists, return it, else
#' create a new one and return it without saving.
#'
#' @return A path to the rank file.
#' @export
#'
#' @examples
fetch_rank_list <- function(ref_condition, treat_condition, out_dir='output'){
  
  rankfile_path <- paste0(c(out_dir, '/', ref_condition, '-', treat_condition, '.rnk'), collapse='')
  
  if (file.exists(rankfile_path)){
    return(read.delim(rankfile_path, sep='\t', row.names = NULL, header=FALSE))
  } else {
    stop("Specified rank file has not been created.")
  }
}



#' Install Genesets from The NCBI
#'
#' Note: not working without login.
#'
#'
#' @param version 
#'
#' @return
#' @export
#'
#' @examples
download_genesets <- function(version='2023.1'){
  get_ncbi_geneset_hyper <- function(gs_name){
    hyper <- paste0(c("https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/", 
                      version, '.Hs/', gs_name, ".all.v", version, ".Hs.symbols.gmt"), collapse='')
    return(hyper)
  }
  
  set_names <- c('h', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'c8')
  
  to_download <- c()
  
  for (set_name in set_names){
    if (! paste0(c(set_name, '.gmt'), collapse='') %in% list.files('genesets')){
      to_download <- c(to_download, set_name)
    }
  }
  
  for (set_name in set_names){
    gs_hyper <- get_ncbi_geneset_hyper(set_name)
    print(gs_hyper)
    download.file(gs_hyper, paste0(c("genesets/", set_name, '.gmt'), collapse=''))
  }
}

# Some useful formatting formatting functions.
s <- function(x){return(as.character(x))}
q <- function(x){return(paste0(c('"', x, '"'), collapse=''))}
twd <- function(x){return(paste0(c(getwd(), '/', x), collapse=''))}

# Helpers to manage storing of past GSEA results.
get_latest_gsea_date <- function(){
  date_to_int <- function(x){
    # Convert the date sting to an integer for comparison 
    # 'YYYY-MM-DD' -> YYYYMMDD
    int_date <- as.numeric(paste0(unlist(strsplit(x, split = '-')), collapse=''))
    return(int_date)
  }
  
  int_to_date <- function(x){
    # Convert an 8 digit date YYYYMMDD -> 'YYYY-MM-DD'
    strdate <- as.character(x)
    strdate <- paste0(c(substr(strdate, 1, 4), substr(strdate, 5, 6),
                        substr(strdate, 7, 8)), collapse = '-')
  }
  
  # Get the latest date string.
  dates <-list.files(twd('gsea_out'))
  int_dates <- unlist(lapply(dates, date_to_int))
  max_int_date <- max(int_dates)
  date <- int_to_date(max_int_date)
  return(date)
}


#' Get the Path to The GSEA Directory Created on The Specified Date
#'
#' @param date 
#' @param safe Whether to check for existence of directory before returning it.
#'
#' @return
#' @export
#'
#' @examples
get_gsea_dir <- function(date='latest', safe=TRUE){
  
  date <- as.character(date)
  
  if (as.character(date) == 'latest'){
    date <- get_latest_gsea_date()
  }
  
  gsea_path <- twd(paste0(c("gsea_out/", as.character(date)), collapse=''))
  if (safe & !file.exists(gsea_path)){
    stop(paste0(c(gsea_path, ' does not exist.'), collapse=''))
  }
  return(gsea_path)
  }


#' Get the Path to the Given Analysis
#'
#' @param reference_cond The reference condition
#' @param compare_cond The condition to be compared.
#' @param date The date on which the comparison was performed.
#'
#' @return A string representation of a path.
#' @export
#'
#' @examples
get_analysis_path <- function(ref_condition, compare_cond, date='latest'){
  
  # Find the matching analysis file name.
  results_dir_contents <- list.files(get_gsea_dir(date))
  analysis_regex <- paste0(c(ref_condition, '_vs_', compare_cond, '.Gsea'), collapse = '')
  results_filename <- results_dir_contents[grepl(analysis_regex, results_dir_contents)]
  
  # Ensure that the correct number of file names were found.
  if (length(results_filename) == 0){
    stop(paste0(c('No matching analysis for', ref_condition, 'and',
                  compare_cond, 'found in', get_gsea_dir(date)), 
                collapse = ' '))
  } else if (length(results_filename) > 1){
    stop(paste0(c('Multiple matching analysis for', ref_condition, 'and',
                  compare_cond, 'found in', get_gsea_dir(date), ':',
                  results_filename), 
                collapse = ' '))
  }
  
  # Return complete path.
  return(paste0(c(get_gsea_dir(date), '/', results_filename), collapse=''))
}


#' Run Gene Set Enrichment Analysis
#'
#' @param rank_path Path to the rank
#' @param seed The random seed with which to do the analysis.
#' @param genesets The sets of genes to check for enrichment.
#'
#' @return
#' @export
#'
#' @examples
get_run_gsea_command <- function(gsea_path, rank_path, seed,
                                 reference_cond, 
                                 compare_cond,
                                 genesets='h.all',
                                 min_set_size=20,
                                 max_set_size=500,
                                 rerun=FALSE,
                                 output_dir=NA){
  
  analysis_name <- paste0(c(reference_cond, '_vs_', compare_cond), collapse = '')
  
  if (is.na(output_dir)){
  output_dir <- get_gsea_dir(date=Sys.Date(), safe=FALSE)
  } else if (output_dir =='gs'){
    output_dir <- get_gsea_dir(date=genesets, safe=FALSE)
  }
  
  
  
  cmd <- paste0(c(
    paste0(c(gsea_path, "/gsea-cli.sh"), collapse=""), "GSEAPreranked",
    
    "-gmx", paste0(c("ftp.broadinstitute.org://pub/gsea/msigdb/human/gene_sets/", s(genesets), ".v2023.1.Hs.symbols.gmt"), collapse=''),
    
    "-collapse",  "Collapse",
    
    "-mode", "Abs_max_of_probes",
    
    "-norm",  "meandiv",
    
    "-nperm", "1000",
    
    "-rnd_seed",  s(seed),
    
    "-rnk", q(twd(s(rank_path))),
    
    "-scoring_scheme",  "weighted",
    
    "-rpt_label", q(analysis_name),
    
    "-chip",  "ftp.broadinstitute.org://pub/gsea/msigdb/human/annotations/Human_Gene_Symbol_with_Remapping_MSigDB.v2023.1.Hs.chip",
    
    "-create_svgs", "false",
    
    "-include_only_symbols","true",
    
    "-make_sets", "true",
    
    "-plot_top_x",  "20",
    
    "-set_max", s(max_set_size),
    
    "-set_min",  s(min_set_size),
    
    "-zip_report", "false",
    
    "-out",  q(output_dir)), collapse=' ')
  
  return(cmd)
}


#' Get The Collapsed Rank File Generated by GSEA
#'
#' @param reference_cond The reference condition.
#' @param compare_cond The comparison condition.
#'
#' @return A data frame containing the ranked list.
#' @export
#'
#' @examples
get_gsea_collapsed_rank_list <- function(ref_condition, compare_cond, 
                                         analysis_name='latest'){
  
  analysis_path <- get_analysis_path(ref_condition, compare_cond, date=analysis_name)
  
  hyphen_title <- paste0(c(ref_condition, '-', compare_cond), collapse = '')
  
  collapsed_rank_path <- file.path(analysis_path, 'edb',
                                   paste0(c(hyphen_title, '_collapsed.rnk'), collapse=''))
  
  return(read.delim(collapsed_rank_path, sep='\t', row.names = NULL, header=FALSE))
}


#' Get the Number of Genes That Were Pruned During the Collapsing Process
#'
#' @param reference_cond The reference condition.
#' @param compare_cond The comparison condition.
#'
#' @return A data frame containing the ranked list.
#' @export
#'
#' @examples
get_gsea_omitted_genes <- function(ref_condition, compare_cond, 
                                   verbose=TRUE, analysis_name='latest'){
  rank_file <- fetch_rank_list(ref_condition, compare_cond, 'output')
  collapsed_rankfile <- get_gsea_collapsed_rank_list(ref_condition, compare_cond,
                                                     analysis_name)
  
  num_missing <- length(rank_file$V1) - length(collapsed_rankfile$V1)
  
  if (verbose){
    print(paste0(c('For', ref_condition, 'vs', compare_cond, 
                   as.character(num_missing),
                   'genes were ommited due to collapsing:'), collapse =' '))
  }
  return(num_missing)
}


#' Get a Path to a File Matching a Pattern in A Directory 
#'
#' @param dir The directory to search.
#' @param patt The pattern to search for.
#'
#' @return
#' @export
#'
#' @examples
get_matching_file <- function(dir, patt){
  all_files <- list.files(dir)
  matching_file <- all_files[grepl(patt, all_files)]
  
  # Ensure that the correct number of file names were found.
  if (length(matching_file) == 0){
    stop(paste0(c('No matching analysis for', patt, 'found in', dir), 
                collapse = ' '))
  } else if (length(matching_file) > 1){
    stop(paste0(c('Multiple matching analysis for', patt, 'found in', dir, ':',
                  matching_file), 
                collapse = ' '))
  }
  return(file.path(dir, matching_file))
}

#' Get a Dataframe Contining the Results From the GSEA Analysis.
#'
#' @param ref_condition 
#' @param compare_cond 
#' @param date 
#'
#' @return
#' @export
#'
#' @examples
get_gsea_results <- function(ref_condition, compare_cond, analysis_name){
  
  analysis_path <- get_analysis_path(ref_condition, compare_cond, analysis_name)
  
  # Extract results by combining positively and negatively enriched genes.
  pos_path <- get_matching_file(analysis_path, 'gsea_report_for_na_pos.*.tsv')
  neg_path <- get_matching_file(analysis_path, 'gsea_report_for_na_neg.*.tsv')
  
  pos_df <- read.delim(pos_path, header=TRUE, sep='\t')
  neg_df <- read.delim(neg_path, header=TRUE, sep='\t')
  
  results <- rbind(pos_df, neg_df)
  
  # Make sure that we have the expected number of gene sets.
  gene_sets_path <- file.path(analysis_path, 'edb', 'gene_sets.gmt')
  expected_num_gene_sets <- R.utils::countLines(gene_sets_path)
  
  if (nrow(results) != expected_num_gene_sets){
    stop(paste0(c(as.character(expected_num_gene_sets - nrow(results)),
                ' gene sets were not analysed by GSEA for unknown reasons.'),
                collapse=''))
  }
  return(results)
}

#' Get a Mapping From Vector A to Vector B
#' 
#' Where a and b are both permutations of each other, with each containing
#' no duplicate entities.
#'
#' @param a 
#' @param b 
#'
#' @return a list of indices, s.t. all(a[inds] == b) == TRUE
#' @export
#'
#' @examples
get_mapping <- function(a, b){
  val_inds <- c()
  for (val in b){
    pos <- which(a == val)
    val_inds <- c(val_inds, pos)
  }
  return(val_inds)
}


get_enrichment_heatmap <- function(results,
                                   n=50,
                                   max_fdr=0.05,
                                   sort_by=c('fdr', 'delta_nes')){
  num_gs <- nrow(results[[1]])
  gs_names <- unlist(results[[1]][1])
  
  # Aggregate NES scores from results into a single matrix.
  nes_matrix <- matrix(ncol=length(results), nrow=num_gs)
  names_matrix <- matrix(ncol=length(results), nrow=num_gs)
  fdr_matrix <- matrix(ncol=length(results), nrow=num_gs)
  row.names(nes_matrix) <- gs_names
  colnames(nes_matrix) <- names(results)
  
  for (i in seq(1, length(results))){
    
    # Use gene set names to map NES and p-vals to the same positions,
    # keep names for confirmation.
    condition_gs_names <- unlist(results[[i]]$NAME)
    condition_nes <- unlist(results[[i]]$NES)
    condition_fdrs <- unlist(results[[i]]$FDR.q.val)
    mapping <- get_mapping(condition_gs_names, gs_names)
    
    nes_matrix[1:num_gs, i] <- condition_nes[mapping]
    names_matrix[1:num_gs, i] <- condition_gs_names[mapping]
    fdr_matrix[1:num_gs, i] <- condition_fdrs[mapping]
    
    if (! all(condition_gs_names[mapping] == gs_names)){
      stop('Aggregating NES scores failed due to mapping error.')
    }
  }
  
  # Filter by FDR threshold.
  min_fdrs <- apply(fdr_matrix, 1, min)
  gt_threshold <- min_fdrs < max_fdr
  nes_matrix <- nes_matrix[gt_threshold,]
  
  # Sort by fdr or delta NES.
  if (any(sort_by == 'fdr')){
    min_fdrs <- min_fdrs[gt_threshold]
    nes_matrix <- nes_matrix[order(min_fdrs, decreasing = FALSE),]
  } 
  
  else if (sort_by == 'delta_nes'){
    
    get_delta <- function(x){
      s <- 0
      for(i in 1:(length(x) - 1)){
        for(j in (i + 1):length(x)){
          s <- s + abs(x[i] - x[j])
        }
      }
      return(s)
    }
    
    deltas <- apply(nes_matrix, 1, get_delta)
    nes_matrix <- nes_matrix[order(deltas, decreasing = TRUE),]
    
  }
  nes_matrix <- nes_matrix[1:min(c(n, nrow(nes_matrix))), ]
  
  heatmap_col = circlize::colorRamp2(c(min(nes_matrix), 0, 
                                       max(nes_matrix)), c("blue", "white", "red"))
  
  
  heatmap <- ComplexHeatmap::Heatmap(nes_matrix,
                          show_row_dend = FALSE, show_column_dend = FALSE, 
                          col=heatmap_col, show_column_names = TRUE, 
                          show_row_names = TRUE,show_heatmap_legend = TRUE, 
                          use_raster=TRUE, row_names_side='left',
                          name='NES', row_names_max_width = unit(15, 'cm'),
                          cluster_columns = FALSE)
  return(heatmap)
}


# Code for getting Hugo Symbols.

# library(biomaRt)
# 
# listMarts()
# ensembl <- useMart("ensembl")
# datasets <- listDatasets(ensembl)
# hg <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
