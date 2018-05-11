#' Retrieve GO Information
#'
#' These function download the GO OBO file and GO annotation file, which are used by the `output_predicted_genes`
#' function.
#'
#' @return NULL, Puts variables GO_GAF and GO in your local environment
#' @export
#'
retrieve_go_information <- function(){
  message("Downloading the GO GAF file, relating genes to GO Terms")
  GO_GAF <- read_tsv("http://geneontology.org/gene-associations/goa_human.gaf.gz",
                     comment = "!", col_names = FALSE) %>%
    suppressWarnings %>%
    select(X3, X5)
  names(GO_GAF) <- c("GENE", "GO_TERM")
  GO_GAF %<>% filter(GO_TERM %in% GO$id) %>% distinct
  GO_GAF <<- GO_GAF
}


#' @export
#' @rdname retrieve_go_information
retrieve_go_OBO <- function(){
  message("Downloading the Gene Ontology (GO) OBO file")
  GO <<- get_ontology("http://purl.obolibrary.org/obo/go.obo")
}

#' Predict disease genes from a seed gene list
#'
#' This function utlizes the structure of the Gene Ontology resource to predict
#' unreported disease genes from a list of seed genes. It accomplishes this task
#' by finding the most granular GO term that a particular gene belongs to, and
#' collecting all other genes belonging to that term. It scores predicted genes
#' based on their frequencies in the related GO terms.
#'
#' By default, only 10 seed genes are used to predict new genes. This can be changed, but increasing
#' this parameter will increase runtime exponentially
#'
#' @param seed_gene_list The tibble outputted by 'output_gene_prioritization'
#' @param n_seeds The number of seed genes to use to extend genes. Defaults to 10
#'
#' @return A tibble
#' @export
#'
#' @examples
#' autism_seed_genes <- output_gene_prioritization("autism")
#' output_predicted_genes(autism_seed_genes)
#'
output_predicted_genes <- function(seed_gene_list, n_seeds = 10) {

  if(sum(names(seed_gene_list) != c("GeneID", "GENE", "reported_SCORE")) > 1)
    stop("The input object must be the output of 'output_gene_prioritization'. Please try again!")

  # Check to make sure the required files are in your local environment
  if(!(exists("GO"))) retrieve_go_OBO()
  if(!(exists("GO_GAF"))) retrieve_go_information()

  seed_genes <- seed_gene_list$GENE %>% .[1:n_seeds]

  score <- seed_gene_list$reported_SCORE
  names(score) <- seed_genes

  # Make sure the gene is in the database before performing extension
  seed_genes %<>% .[. %in% GO_GAF$GENE]

  # For each gene in the seed gene list, find related genes given the following protocol
  extended_genes_from_seeds_list <- lapply(seed_genes, function(gene) {

    # Get the number of children for each matched GO term
    n_children <-
      GO$children[filter(GO_GAF, GENE == gene)$GO_TERM] %>% map_dbl(length)

    # Get the most granular GO terms that the gene belongs to. In most cases, this is a GO term
    # with no children
    granular_GOs <- n_children[n_children == min(n_children)] %>% names

    # For each granular GO term, collect a list of other genes that are members. Assign those members
    # a score equal to the seed gene's score, divided by the square root of the number of genes in the GO term
    # This way, GO terms which contain a huge number of genes won't be as informative
    extended_granular_GOs <- lapply(granular_GOs, function(GO) {
      temp <- GO_GAF %>% filter(GO_TERM == GO)
      temp %>% mutate(SCORE = score[gene] / sqrt(nrow(temp))) %>% distinct
    }
      )

    # Aggregate genes by summing their scores across all GO terms. If a gene is present in many of the GO
    # terms that the seed gene is in, it will be given a higher score.
    extended_genes <- extended_granular_GOs %>%
      do.call("rbind", .) %>%
      select(GENE, SCORE) %>%
      group_by(GENE) %>%
      summarise(SCORE = sum(SCORE)) %>%
      arrange(desc(SCORE))
  })

  # In the list of genes produced from each seed gene, combine them by summing their scores
  # Normalize the scores the same way as done previously, by dividing through by the max
  extended_genes_from_seeds_list %>%
    do.call("rbind", .) %>%
    group_by(GENE) %>%
    summarise(SCORE = sum(SCORE)) %>%
    arrange(desc(SCORE)) %>%
    mutate(SOURCE = ifelse(GENE %in% seed_genes, "seed_gene", "predicted_gene")) %>%
    mutate(SCORE = SCORE / SCORE[1])


}
