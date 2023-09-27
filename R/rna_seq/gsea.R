

# === GSEA Directory Management Funcitons ===

#' Get the Path to the Given Analysis
#'
#' @param reference_cond The reference condition
#' @param compare_cond The condition to be compared.
#' @param genesets The genesets with which the analysis was conducted.
#' @param parent_dir The directory containing the analysis sub directories.
#' 
#' @return A string representation of a path.
#' @export
#'
#' @examples
get_analysis_path <- function(ref_condition, compare_cond, genesets, 
                              parent_dir='gsea_out'){
  containing_path <- file.path(parent_dir, genesets)
  
  # Find the matching analysis file name.
  results_dir_contents <- list.files(containing_path)
  analysis_regex <- paste0(c(ref_condition, '_vs_', compare_cond, '_', genesets, '.Gsea'), collapse = '')
  results_filename <- results_dir_contents[grepl(analysis_regex, results_dir_contents)]
  print(results_filename)
  
  # Ensure that the correct number of file names were found.
  if (length(results_filename) == 0){
    stop(paste0(c('No matching analysis for', ref_condition, 'and',
                  compare_cond, 'found in', containing_path), 
                collapse = ' '))
  } else if (length(results_filename) > 1){
    stop(paste0(c('Multiple matching analysis for', ref_condition, 'and',
                  compare_cond, 'found in', containing_path, ':',
                  results_filename), 
                collapse = ' '))
  }
  
  # Return complete path.
  return(paste0(c(containing_path, '/', results_filename), collapse=''))
}

#' Get Whether a Given GSEA Directory Had Been Created
#'
#' @param reference_cond The reference condition
#' @param compare_cond The condition to be compared.
#' @param genesets The genesets with which the analysis was conducted.
#' @param parent_dir The directory containing the analysis sub directories.
#'
#' @return A string representation of a path.
#' @export
#'
#' @examples
analysis_exists <- function(ref_condition, compare_cond, geneset, 
                            parent_dir='gsea_out'){
  containing_path <- file.path(parent_dir, geneset)
  
  # Find the matching analysis file name.
  results_dir_contents <- list.files(containing_path)
  analysis_regex <- paste0(c(ref_condition, '_vs_', compare_cond, '_', geneset, '.Gsea'), collapse = '')
  results_filename <- results_dir_contents[grepl(analysis_regex, results_dir_contents)]
  print(results_filename)
  
  # Ensure that the correct number of file names were found.
  return(length(results_filename) != 0)
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
#' @param reference_cond The reference condition
#' @param compare_cond The condition to be compared.
#' @param genesets The genesets with which the analysis was conducted.
#' @param parent_dir The directory containing the analysis sub directories.
#'
#' @return
#' @export
#'
#' @examples
get_gsea_results <- function(ref_condition, compare_cond, genesets){
  
  analysis_path <- get_analysis_path(ref_condition, compare_cond, genesets)
  
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


#' Get all GSEA Results from the Given Directory
#' 
#' @param gsea_out_dir The output directory containing all of the GSEA results.
#' This directory must only contain directories containing GSEA results.
#' 
#'  |gsea_out_dir
#'  |--- Genesets 1
#'  |  |--- GSEA Result from Condition 1
#'  |  |--- GSEA Result from Condition 2
#'  |--- Genesets 2
#'  |  |--- GSEA Result from Condition 1
#'  |  |--- GSEA Result from Condition 2
#'  |--- h.all
#'  |  |--- D81WT_vs_D81A_h.all
#'  |  |--- D81WT_vs_D81N_h.all
#' 
#' @return A list containing a list of result dataframes.
#' 
#' === Example <gsea_result> Structure ===
#' 
#'|--- Genesets 1
#'|  |--- GSEA Result from Condition 1
#'|  |--- GSEA Result from Condition 2
#'|--- Genesets 2
#'|  |--- GSEA Result from Condition 1
#'|  |--- GSEA Result from Condition 2
#'|--- h.all
#'|  |--- D81WT_vs_D81A_h.all
#'|  |--- D81WT_vs_D81N_h.all
#' 
#' Where each GSEA result contains the values reported in the for both of
#' gsea_report_for_neg_.tsv and gsea_report_for_pos_.tsv.
#' 
#' The <reference>_<comparison>_<genesets> name order is important for 
#' downstream interpretation.
#' 
#' 
get_all_gsea_results <- function(gsea_out_dir='gsea_out', verbose=TRUE){
  
  all_gsea_results <- list()
  for (gs_dir in list.files(gsea_out_dir)){
    
    print(paste0(c('|--\t', gs_dir), collapse=''))
    all_gsea_results[[gs_dir]] <- list()
    
    res_path <- file.path(gsea_out_dir, gs_dir)
    for (gsea_res in list.files(res_path)){
      print(gsea_res)
      
      components <- partition_comparison_name(gsea_res)
      
      all_gsea_results[[gs_dir]][[gsea_res]] <- get_gsea_results(components[1],
                                                     components[2],
                                                     components[3])
      
      analysis_name <- paste0(c(components[1], '_vs_', components[2], '_',
                                components[3]), collapse = '')
      print(paste0(c('|\t|--\t', analysis_name), collapse=''))
    }
  }
  return(all_gsea_results)
}

#' Get the Components of a GSEA Comparison from its Formatted name
#' 
#' @param gsea_result_name character formatted 
#' <ref_cond>_vs_<compare_cond>_<geneset>.Gsea...
#' 
#' @return Length 3 named character vector with components equal to  <ref_cond>, 
#' <compare_cond>, and <geneset>.
#' 
partition_comparison_name <- function(gsea_result_name){
  s1 <- unlist(strsplit(gsea_result_name, '.Gsea'))
  s2 <- unlist(strsplit(s1[1], '_'))
  genesets <- s2[length(s2)]
  
  s3 <- paste0(s2[-length(s2)], collapse='_')
  s4 <- unlist(strsplit(s3, '_vs_'))
  ref_cond <- s4[1]
  compare_cond <- s4[2]
  
  rtrn <- c(ref_cond, compare_cond, genesets)
  names(rtrn) <- c('reference_cond', 'compare_cond', 'genesets')
  return(rtrn)
}


# === Functional GSEA Analysis Code ===

#' Generate a .rnk File Using The Given DE List
#'
#' @param de_list List of DE genes with columns $gene_name, $logFC, and
#'  $adj.P.Val or $P.Val.
#' @param ref_condition The name of the reference (often control) condition.
#' @param treat_codition The name of the treatment condition.
#' @param out_dir The directory in which to place the ranked list.
#' @param return_dataframe If the specified rank file exists, return it, else
#' create a new one and return it without saving.
#' @param mult_hyp_testing Whether to extract multiple hypothesis testing P values,
#' as opposed to regular p values.
#'
#' @return A path to the rank file.
#' @export
#'
#' @examples
get_rank_list <- function(de_list, ref_condition, treat_condition, out_dir='output',
                          mult_hyp_testing=TRUE){
  
  rankfile_path <- paste0(c(out_dir, '/', ref_condition, '-', treat_condition, '.rnk'), collapse='')
  
  rnk <- data.frame(gene_name=de_list$gene_name)
  
  dir_pvals <- numeric()
  for (i in seq_along(de_list$adj.P.Val)){
    logfc <- de_list$logFC[i]
    
    if (mult_hyp_testing){
      pval <- de_list$adj.P.Val[i]
    } else {
      pval <- de_list$P.Value[i]
    }
    
    dir_pvals[i] <- -log10(pval) * sign(logfc)
  }
  
  de_list$rank <- dir_pvals
  
  rnk$dir.adj.P.Val <- dir_pvals
  
  rnk <- rnk[order(rnk$dir.adj.P.Val, decreasing=TRUE),]
  
  write.table(rnk, rankfile_path, sep='\t', row.names=FALSE, col.names = FALSE)
  
  return(rankfile_path)
}



# Some useful formatting formatting functions.
s <- function(x){return(as.character(x))}
q <- function(x){return(paste0(c('"', x, '"'), collapse=''))}
twd <- function(x){return(paste0(c(getwd(), '/', x), collapse=''))}



#' Run Gene Set Enrichment Analysis
#' 
#' @param gsea_path Path to the folder containing command line install of GSEA.
#' 
#' @param rank_path Path to rankfile.
#' 
#' @param seed The seed to use to generate the analyses, for reproducability.
#' 
#' @param reference_cond The reference/control condition.
#' 
#' @param compare_cond The treatment/comparison condition. A positive enrichment
#'  score indicates a positive enrichment in this condition.
#'  
#' @param genesets The genesets to use, as referenced by the Broad MSIGDB.
#' 
#' @param min_set_size The minimum size of any given gene set to be analysed.
#' 
#' @param max_set_size The maximum size of any given gene set to be included.
#' 
#' @param rerun Whether to re-run the analysis if a matching analysis already 
#'  exists at a given path. Iff false, creates a second directory for the new
#'  analysis with a slightly different name.
#'  
#' @param output_dir The output directory. Results will be placed in 
#'  <output_dir>/<genesets>/<ref_condition>_vs_<compare_cond>_<genesets>.Gsea
#'
#' @param rank_path Path to the rank
#' @param seed The random seed with which to do the analysis.
#' @param genesets The sets of genes to check for enrichment.
#'
#' @return
#' @export
#'
#' @examples
#' 
#' get_run_gsea_command("/n/data1/dfci/pedonc/kadoch/ben/apps/GSEA_4.3.2", 
#'                       "output/dummy_rnk.rnk", 52, 'DMSO', 'FHD286', 
#'                       genesets='h.all')
#' 
#' 
#' 
get_run_gsea_command <- function(gsea_path, 
                                 rank_path, 
                                 seed,
                                 reference_cond, 
                                 compare_cond,
                                 genesets='h.all',
                                 min_set_size=20,
                                 max_set_size=500,
                                 rerun=FALSE,
                                 output_dir='gsea_out',
                                 run=FALSE){
  
  analysis_name <- paste0(c(reference_cond, '_vs_', compare_cond, '_',
                            genesets), collapse = '')
  
  if (! rerun & analysis_exists(reference_cond, compare_cond, genesets,
                                parent_dir=output_dir)){
    warning(paste0(c(analysis_name, 'already exists in', output_dir, '/',
                     'genesets. Avoiding re-run.'), collapse=' '))
    return('echo Avoiding Re-run')
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
    
    "-create_svgs", "true",
    
    "-include_only_symbols","true",
    
    "-make_sets", "true",
    
    "-plot_top_x",  "20",
    
    "-set_max", s(max_set_size),
    
    "-set_min",  s(min_set_size),
    
    "-zip_report", "false",
    
    "-out",  q(c(output_dir, '/', genesets))), collapse=' ')
  
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
get_gsea_collapsed_rank_list <- function(ref_condition, compare_cond, genesets){
  
  analysis_path <- get_analysis_path(ref_condition, compare_cond, genesets)
  
  print(analysis_path)
  
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
get_gsea_omitted_genes <- function(ref_condition, compare_cond, genesets,
                                   rank_file_parent='output', 
                                   gsea_parent='gsea_out', verbose=TRUE){
  rank_file <- fetch_rank_list(ref_condition, compare_cond, rank_file_parent)
  collapsed_rankfile <- get_gsea_collapsed_rank_list(ref_condition, compare_cond,
                                                     analysis_name, gsea_parent)
  
  num_missing <- length(rank_file$V1) - length(collapsed_rankfile$V1)
  
  if (verbose){
    print(paste0(c('For', ref_condition, 'vs', compare_cond, 
                   as.character(num_missing),
                   'genes were ommited due to collapsing:'), collapse =' '))
  }
  return(num_missing)
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


# === Figure Generation ===


#' Produce a Heatmap Visualising the Gene Set Enrichment REsults
#'
#' @param results A list of GSEA results, where each entry corresponds to the
#' concatentation of gsea_report_for_na_pos.tsv and 
#' gsea_report_for_na_neg.tsv generated from a single GSEA analysis.
#' 
#' @param n The maximum number of gene sets to display.
#' 
#' @param max_fdr In order for a gene set to be included as a row in the matrix,
#' it must be enriched with it's fdr < max_fdr in at least oen of the 
#' conditions.
#' 
#' @param sort_by One of "fdr" or "delta_nes".
#'  Will include the gene sets with the top abs(NES) or FDR.
#'
#' @return A heatmap object.
#' @export
#'
#' @examples
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





#' Generate a Bubble-Lattice Plot Displaying GSEA Results
#'
#' @param all_results A list of lists of the structure described as in the 
#'  function *get_all_gsea_results*.
#'  
#' @param genesets_of_interest A list of geneset names, to be displayed in
#'  order on the y axis of the 
#' 
#' @param comparisons_of_interest The names of the comparisons of interest.
#'  A list of length 2 vectors, with the first entry being the reference
#'  condition, and the second entry being the comparison condition.
#'
#' @return A bubble-lattice plot.
#' @export
#'
#' @examples
bubble_lattice_plot <- function(all_results, 
                                genesets_of_interest, 
                                comparisons_of_interest,
                                fdr_threshold){
  num_rows <- length(genesets_of_interest) * length(lines)
  
  formatted_results <- data.frame(DataSet=c(),
                                  GeneSet=c(),
                                  NES=numeric(),
                                  isSig=logical())
  
  # Iterate throughout all_results, looking for entries of one of 
  # <genesets_of_interest> for any of <comparisons_of_interest>
  for (genesets in names(all_results)){
    comparison_results <- all_results[[genesets]]
    
    for (comparison in names(comparison_results)){
      parts <- partition_comparison_name(comparison)
      reference <- parts[['reference_cond']]
      compare <- parts[['compare_cond']]
      
      # This isn't a comparison we care about, so skip it.
      comp <- list(c(reference, compare))
      if (! comp %in% comparisons_of_interest){
        next
      }
      
      comp_data_frame <- comparison_results[[comparison]]
      comparison <- paste0(c(reference, '_vs_', compare), collapse='')
      
      for (row in 1:nrow(comp_data_frame)){
        geneset <- comp_data_frame$NAME[row]
        
        if (!geneset %in% genesets_of_interest){
          next
        }
        
        nes <- comp_data_frame$NES[row]
        fdr <- comp_data_frame$FDR.q.val[row]
        is_sig <- fdr < fdr_threshold
        
        new_row <- c(comparison, geneset, nes, is_sig)
        print(new_row)
        
        formatted_results <- rbind(formatted_results, new_row)
      }
    }
  }
  
  colnames(formatted_results) <- c("DataSet", "GeneSet", "NES", "isSig")
  
  formatted_results$NES <- as.numeric(formatted_results$NES)
  
  library(colorspace)
  library(RColorBrewer)
  my_color_scale <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
  
  data <- formatted_results
  # Build the ggplot object
  plt <- ggplot(data, aes(x=DataSet, y=GeneSet, size=abs(NES))) +
    geom_point(aes(fill=NES, alpha=isSig, color=isSig), pch=21, stroke = 1) + 
    theme_light() +
    theme(axis.text.x = element_text(angle = 45, hjust=1)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_fill_gradientn(colors = my_color_scale(100), limits = c(-4,4)) +
    scale_color_manual(values = c("grey", "black"), name=paste0(c("FDR < ", as.character(fdr_threshold)), collapse = '')) +
    scale_alpha_manual(values = c(0.3, 1)) +
    scale_size_continuous(trans = "log10", range = c(0,8))
  return(plt)
}




# Code for getting Hugo Symbols.

# library(biomaRt)
# 
# listMarts()
# ensembl <- useMart("ensembl")
# datasets <- listDatasets(ensembl)
# hg <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
