retrieve_go_information <- function(){
  GO_GAF <- read_tsv("http://geneontology.org/gene-associations/goa_human.gaf.gz",
                     comment = "!", col_names = FALSE) %>%
    suppressWarnings %>%
    select(X3, X5)
  names(GO_GAF) <- c("GENE", "GO_TERM")
  GO_GAF %<>% filter(GO_TERM %in% GO$id)
  GO_GAF <<- GO_GAF
}

predict_genes_from_seed_list <- function(seed_gene_list) {

  if(!(exists("GO_GAF") & exists("CTD"))) retrieve_go_information()
  seed_genes <- seed_gene_list$GENE

  score <- seed_gene_list$reported_SCORE
  names(score) <- seed_genes

  # Make sure the gene is in the database before performing extension
  seed_genes %<>% .[. %in% GO_GAF$GENE]

  # For each gene in the seed gene list, find related genes given the following protocol
  extended_genes_from_seeds_list <- lapply(seed_genes, function(gene) {
    n_children <-
      GO$children[filter(GO_GAF, GENE == gene)$GO_TERM] %>% map(length) %>% unlist

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
