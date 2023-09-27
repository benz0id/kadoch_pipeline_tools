#' Get a Dataframe Mapping a Given Site to its Nearest Gene.
#'
#' @param nearest_gene_map_path Path to a bed file annotated with the nearest
#'  gene.
#'
#' @return A dataframe with two columns, where each row corresponds to a site
#'  and its nearest gene.
#' @export
#' 
#' Note that the returned dataframe may be missing some sites, which were lost
#' when trying to find the nearest gene.
#'
#' @examples
get_site_to_gene_map <- function(nearest_gene_map_path){
  
  nearest_gene <- read.csv2(nearest_gene_map_path, sep='\t')
  sites <- c()
  for (row in seq(1, nrow(nearest_gene))){
    site <- paste0(c(nearest_gene$chr[row], ':', 
                     as.character(nearest_gene$start[row]), '-',
                     as.character(nearest_gene$stop[row])), collapse='')
    sites <- c(sites, site)
  }
  sites_to_nearest_gene <- data.frame(site=sites, gene=nearest_gene$Gene)
  return(sites_to_nearest_gene)
}



#' Find the Nearest Gene for a Given Site
#'
#' @param site A given genomic location, formatted <contig>:<start>-<stop>.
#' @param site_to_gene_map A dataframe with columns "site" and "gene", mapping
#'  a given site to its nearest gene.
#' @param verbose 
#'
#' @return The nearest gene to the given site. Character vector of length 1
#' @export
#'
#' @examples
site_to_gene_name <- function(site, site_to_gene_map, verbose=FALSE){
  site_pos <- which(site_to_gene_map$site == site)
  
  if (length(site_pos) != 1){
    if (verbose){
      print('no gene found')
      print(site)
      print(sum(site_pos))
    }
    return(NA)
  }
  
  site_gene <- site_to_gene_map$gene[site_pos]
  
  return(site_gene)
}

#' Get all of the Nearest Genes for all of The Given Sites
#'
#' @param sites A vector of site strings, <contig>:<start>-<stop>.
#' @param site_to_gene_map A dataframe mapping a given site to its nearest genes.
#' @param verbose 
#'
#' @return Character vector of length equal to that of sites. Sites for which no
#'  gene is found will have an NA placed at the corresponding position.
#' @export
#'
#' @examples
get_nearest_genes <- function(sites, site_to_gene_map, verbose=FALSE){
  nearest_genes <- sapply(sites, site_to_gene_name, site_to_gene_map=site_to_gene_map)
  if (verbose){
    missed_sites <- sites[is.na(nearest_genes)]
    print(paste0(c('Sites missed: ', missed_sites), collapse='\n\t'))
  }
  return(nearest_genes)
}



#' Convert a Dataframe of DE Sites into a dataframe of DE Nearest Genes.
#' 
#' Given a dataframe containing columns "gene_names" - which contains sites 
#' formatted <contig>:<start>-<stop> and a column "adj.P.Val", replaces all of
#' the given sites with their nearest genes.
#' 
#' Rows with duplicate gene names are also removed, with preference being 
#' given to rows with lower adj.P.Val.
#' 
#' I know it's weird that a column that is called
#' gene_names contains sites I'm sorry I'll fix this eventually.
#'
#' @param de_list A dataframe containing columns "gene_name" and "adj.P.Val",
#'  containing sites and adjusted P value, respective.
#' @param site_to_gene_map A dataframe mapping sites to nearest genes.
#'
#' @return A copy of <de_list>, with the sites in "gene_name" replaced by actual
#'  gene names. 
#'  
#'    Sites without entries in <site_to_gene_map> are pruned.
#'  
#'    If gene names are duplicated, only the row with the greatest
#'  p_value will be kept. 
#'  
#' @export
#'
#' @examples
convert_de_to_gene_names <- function(de_list, site_to_gene_map, verbose=TRUE){
  gene_names <- get_nearest_genes(de_list$gene_name, site_to_gene_map)
  de_list$sites <- de_list$gene_name
  de_list$gene_name <- gene_names
  
  nas <- is.na(de_list$gene_name)
  de_list <- de_list[!nas, ]
  
  de_list <- de_list[order(de_list$adj.P.Val), ]
  
  dups <- duplicated(de_list$gene_name)
  de_list <- de_list[!dups,]
  
  print(paste0(c(as.character(sum(nas)), 'na genes pruned and', as.character(sum(dups)),
                 'duplicates removed.'), collapse=' '))
  return(de_list)
}













