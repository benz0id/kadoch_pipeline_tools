#' Get a Dataframe Mapping Any* Given ENSEMBL ID to it's Correspoding HUGO Symbol
#'
#' @return
#' @export
#'
#' @examples
get_ensembl_map <- function(){
  mart <- biomaRt::useEnsembl("ensembl", "hsapiens_gene_ensembl")
  symb <- biomaRt::keys(org.Hs.eg.db::org.Hs.eg.db, "SYMBOL")
  map <- biomaRt::getBM(c("ensembl_gene_id","hgnc_symbol"), "hgnc_symbol", symb, mart)
  return(map)
}

get_ensembl_id <- function(gene_symbol, map, verbose=TRUE){
  matching <- map$hgnc_symbol == gene_symbol
  
  if (verbose & all(!matching)){
    print(paste0(c("Failed to find any HUGO gene symbols matching", 
                   gene_symbol), collapse = ' '))
    return(NA)
  } else if (sum(matching) > 1){
    print(paste0(c("Found multiple HUGO gene symbols matching", 
                   gene_symbol, ':', map$ensembl_gene_id[matching]), 
                 collapse = ' '))
    return(map$ensembl_gene_id[matching][1])
  }
  
  return(map$ensembl_gene_id[matching])
}

#' Get Ensembl IDs for Each of the Given Gene Symbols
#'
#' @param gene_symbols 
#' @param ensembl_to_hugo 
#' @param verbose 
#'
#' @return 
#' @export
#'
#' @examples
symbols_to_ensembl <- function(gene_symbols, ensembl_to_hugo, verbose=TRUE,
                               update_symbols=TRUE){
  if (update_symbols){
    gene_symbols <- limma::alias2Symbol(gene_symbols, species = 'Hs',
                                        expand.symbols = FALSE)
  }
  ensembl_ids <- sapply(gene_symbols, get_ensembl_id, map=ensembl_to_hugo, verbose=verbose)
  return(ensembl_ids)
}


#' Get an Updated List of Gene Symbols
#' 
#' Use limma to search NCBI databases for past aliases of gene names. If a gene
#' cannot be confidently assigned to a new symbol, keep the old one.
#'
#' @param gene_symbols The gene symbols to be updated.
#' @param output_folder The folder to output the results to.
#' @param save_omitted_genes Whether to save the genes that were ommited to the 
#' output folder.
#'
#' @return A list of updated gene names.
#' @export
#'
#' @examples
update_symbols <- function(gene_symbols, output_folder='output', 
                           save_omitted_genes=TRUE, verbose=TRUE){
  new_symbols <- limma::alias2SymbolTable(gene_symbols)
  not_updateable <- is.na(new_symbols)
  
  updated_symbols <- gene_symbols
  updated_symbols[!not_updateable] <- new_symbols[!not_updateable]
  
  if (save_omitted_genes){
    write.csv(data.frame(updated_symbols[not_updateable]), 
              paste0(c(output_folder, '/unupdateable_genes.csv'), collapse=''))
  }
  if (verbose){
    print(paste0(c(as.character(sum(updated_symbols != gene_symbols)), '/', 
            as.character(length(gene_symbols)), ' gene symbols were updated, ',
            as.character(sum(not_updateable)), 
            ' failed to map to a gene symbol.'), collapse=''))
  }
  
  return(updated_symbols)
}


#' Remove NA and Duplicated Rows From a Dataframe
#'
#' @param data The dataframe to remove rows from.
#' @param attribute_vector The attribute to filter by.
#' @param names_vector The name of each row. Only required if verbosity is enabled.
#' @param remove_na Whether to remove NA values. 
#' @param remove_dups Whehter to remove duplicated values
#' @param verbose 
#'
#' @return
#' @export
#'
#' @examples
filter <- function(data, attribute_vector, names_vector=NA, remove_na=TRUE, 
                   remove_dups=TRUE, verbose=TRUE){
  
  if (verbose & all(is.na(names_vector))){
    stop("Names vector must be present if verbosity requested.")
  }
  
  to_remove <- logical(length(attribute_vector))
  
  if (remove_na){
    to_remove <- to_remove | is.na(attribute_vector)
    if (verbose){
      num_na <- sum(is.na(attribute_vector))
      print(paste0(c(as.character(num_na), ' NAs found.'),
                   collapse= ''))
    }
  }
  
  if (remove_na){
    to_remove <- to_remove | duplicated(attribute_vector)
    if (verbose){
      num_dups <- length(attribute_vector) - length(unique(attribute_vector))
      print(paste0(c(as.character(num_dups), ' duplicates found.'),
                   collapse= ''))
    }
  }
  
  if (verbose & any(to_remove)){
    print(c('Removed:', names_vector[to_remove]))
  } else if (verbose){
    print(c('No entries removed.'))
  }
  
  return(data[!to_remove,])
  
}