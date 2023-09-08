source("dge_helpers.R")


get_mean_expression <- function(norm_cpms, treatment_name, groups, 
                                col_offset=1){
  
  # Extract treatment columns and return their average
  treatment_cols <- which(groups == treatment_name) + col_offset
  treatment_df <- norm_cpms[treatment_cols]
  
  # Convert from log2 cpm to cpm.
  treatment_cpm <- apply(treatment_df, MARGIN=2, FUN=function(x){
    return(2 ^ x)
  })
  averages <- apply(treatment_cpm, MARGIN=1, FUN=mean)
  avgs_df <- data.frame(gene_name=norm_cpms$gene_name, 
                        average_reads=averages)
  
  return(avgs_df)
}

#' Get A Matrix Scoring the Synergy between genes.
#'
#' @param norm_cpms The normalized counts per million.
#' 
#' @param control_name The name of the control condition.
#' 
#' @param t1_name The name of the first treatment condition.
#' 
#' @param t2_name The name of the second treatment condition.
#' 
#' @param t1_and_t2_name The name of the joint treatment condition.
#' 
#' @param groups Vector of treatment conditions, one for each column in 
#' <norm_cpms>.
#' 
#'
#' @return
#' @export
#'
#' @examples
get_synergy_matrix <- function(norm_cpms, control_name, t1_name, t2_name, 
                               t1_and_t2_name, groups, num_desc_columns=1){
  
  # Collect average gene expression for each condition.
  ave_expression <- data.frame(gene_name=norm_cpms$gene_name)
  for (group_name in c(control_name, t1_name, t2_name, t1_and_t2_name)){
    mean_expr_df <- get_mean_expression(norm_cpms, group_name, groups)
    if (! all(mean_expr_df$gene_name == norm_cpms$gene_name)){
      stop("Detected misalignemnt between gene names.")
    }
    ave_expression[[`group_name`]] <-  mean_expr_df$average_reads
  }
  
  # Search for non additive iterations.
  ave_expr_matrix <- as.matrix(ave_expression[1:nrow(ave_expression),
                                              2:5])
  rownames(ave_expr_matrix) <- ave_expression$gene_name
  
  # Extract Expression under each condition.
  cond_expr <- ave_expr_matrix[1:nrow(ave_expr_matrix),
                               1:4]
  colnames(cond_expr) <- c('C', 'A', 'B', 'AB')
  
  # Compute synergy.
  fc_cond_expr <- cbind(cond_expr,
                        'A/C'=cond_expr[,'A'] / cond_expr[,'C'],
                        'B/C'=cond_expr[,'B'] / cond_expr[,'C'],
                        'AB/C'=cond_expr[,'AB'] / cond_expr[,'C'])
  
  synergy_matrix <- cbind(fc_cond_expr, 
                        'S'=fc_cond_expr[, 'AB/C'] 
                        - fc_cond_expr[, 'A/C'] 
                        - fc_cond_expr[, 'B/C'] + 1)
  
  colnames(synergy_matrix)[1:4] <- c(control_name, t1_name, t2_name, t1_and_t2_name)
  
  return(synergy_matrix)
                          
}


#' Plot Synergy Matrices
#'
#' @param synergy_matrix The synergy matrix containing the synergies of interest.
#' 
#' @param n The number of genes to visualise on the barplot.
#' 
#' @param do_top Whether to present the genes with the highest (or lowest)
#' synergy.
#' 
#' @param do_log Whether to transform the fold change values into log fold change
#' values.
#'
#' @return The plot
#' @export
#'
#' @examples
synergy_barplot <- function(synergy_matrix, n, do_top=TRUE, do_log=TRUE){
  
  # Subset synergy matrix.
  synergy_matrix <- synergy_matrix[order(synergy_matrix[, 'S'], 
                                              decreasing=do_top),][1:n,]
  
  if (do_log){
    fun <- log2
  } else {
    fun <- function(x){return(x)}
  }
  
  control_name <- colnames(synergy_matrix)[1]
  
  data <- data.frame(test=c(1, 2, 3))
  for (i in seq(n, 1, -1)){
    gene_name <- rownames(synergy_matrix)[i]
    data[[`gene_name`]] <- c(fun(synergy_matrix[i,'A/C']), 
                             fun(synergy_matrix[i,'B/C']), 
                             fun(synergy_matrix[i,'AB/C']))
  }
  data <- data[2:ncol(data)]
  rownames(data) <- colnames(synergy_matrix)[2:4]
  
  if (do_log){
    ylab <- paste0(c("Log2 Fold Change Over ", control_name), collapse = '')
  } else {
    ylab <- paste0(c("Fold Change Over ", control_name), collapse = '')
  }
  
  if (! do_top){
    title <- 'Changes in Expression Amoung Groups With the Lowest Synergy'
  } else {
    title <- 'Changes in Expression Amoung Groups With the Highest Synergy'
  }
  
  colours <- c("lightblue", "pink", "plum1")
  
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  
  plt <- barplot(as.matrix(data),
          main=title,
          
          # setting y label only
          # because x-label will be our
          # barplots name
          ylab=ylab,
          
          # to plot the bars vertically
          beside=TRUE,
          col=colours,
          cex.names=0.7
  )
  
  legend("topright", 
         legend = rownames(data), 
         fill = colours,
         inset=c(-0.3,0),
         cex=0.7)
  
  return(plt)
}


additive_dotplot <- function(synergy_matrix, n, do_top=TRUE){
  
  # Subset synergy matrix.
  synergy_ranked <- synergy_matrix[order(synergy_matrix[, 'S'], 
                                         decreasing=do_top),][1:n,]
  synergy_ranked <- cbind(synergy_ranked, 1:n)
  colnames(synergy_ranked)[ncol(synergy_ranked)] <- 'rank'
  synergy_ranked <- as.data.frame(synergy_ranked)
  
  if (do_top){
    pref <- 'Most'
  } else {
    pref <- 'Least'
  }
  
  title <- paste0(c('Synergy Amoung ', pref, ' Synergistic Genes Ordered by Rank'), collapse = '')
  
  
  ggplot(synergy_ranked, aes(rank, S)) + 
    geom_point() + 
    ylab("Synergy (Fold Change)") +
    xlab('Rank') +
    labs(title=title)
    
}





















