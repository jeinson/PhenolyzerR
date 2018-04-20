#' This function uses the public ontology_search function to try both databases
#' @export
disease_extension <- function(term, exact_match = FALSE) {

  # If the ontology files aren't in the global environment, get them
  if(!(exists("doid") & exists("ctd"))) retrieve_ontologies()

  term <- term

  doid_results <- ontology_search(term, doid)

  message("Finished searching the Disease Ontology Database")
  if (length(doid_results > 0)) {
    doid_results <-
      cbind(doid_results, rep("DISEASE_ONTOLOGY", length(doid_results)))
  }
  else
    doid_results <-
    tibble(a = character(0), b = character(0)) # If no results found, make a dummy variable


  ctd_results <- ontology_search(term, ctd)

  message("Finished searching the CTD Medic Database")
  if (length(ctd_results) > 0) {
    ctd_results <-
      cbind(ctd_results, rep("CTD_DISEASE", length(ctd_results)))
  }
  else
    ctd_results <-
    tibble(a = character(0), b = character(0))  # If no results found, make a dummy variable

  output <- as.tibble(rbind(doid_results, ctd_results))
  names(output) <- c("Disease", "Source")
  output$Source <- output$Source %>% as.factor
  output

}

#' @export
output_gene_prioritization <- function(disease_term) {
  useless_words <- c("disease", "syndrome")
  # Split the query string into individual words and get rid of common "useless words"
  term_k <-
    disease_term %>% str_split(" ") %>% unlist %>% discard( ~ .x %in% useless_words)

  if (length(term_k) > 1){print(paste("Splitting your query into", length(term_k), "searches"))}

  # Make sure databases are all there
  if(!(exists("doid") & exists("ctd"))) build_gene_disease_reference()

  # Run extension on each Term_k
  diseases_k <- lapply(term_k, disease_extension)
  if (max(sapply(diseases_k, nrow)) == 0)
    stop("No matching disease terms found. Please check your spelling!!")

  # Clean each Term_k and remove duplicates
  for (i in 1:length(diseases_k)) {
    diseases_k[[i]]$Disease %<>% cleaner
  }
  diseases_k <-
    lapply(diseases_k, function(x)
      x %>% mutate(Disease = Disease %>% tolower %>% trimws) %>% distinct)
  names(diseases_k) <- term_k

  # Make sure DB_COMPILED_DISEASE_SCORES is in the global environnment
  if(!exists("DB_COMPILED_GENE_DISEASE_SCORES"))

  # Match each disease in diseases_k with the precompiled database
  disease_scores <- lapply(diseases_k, function(x) {
    DB_COMPILED_GENE_DISEASE_SCORES %>% filter(tolower(.$DISEASE) %>% map_lgl( ~ . %in% x$Disease))
  })

  if (max(sapply(disease_scores, nrow)) == 0)
    stop("No disease / gene relationships found. Please adjust your query!!")

  # Aggregate the disease_k list
  out <-
    lapply(disease_scores, function(disease_k)
    {
      disease_k %>% group_by(GENE) %>% summarise(total_SCORE = sum(SCORE) / nrow(disease_k))
    }) %>%
    do.call("rbind", .) %>%
    group_by(GENE) %>%
    summarise(reported_SCORE = sum(total_SCORE)) %>%
    arrange(desc(reported_SCORE))

  out$reported_SCORE <- out$reported_SCORE / out$reported_SCORE[1]
  out <- left_join(out, HUMAN_GENE_ID %>% select(GeneID, Symbol), by = c("GENE" = "Symbol"))
  out %<>% select("GeneID", "GENE", "reported_SCORE")

  out
  # Possible place to put in a wordcloud function, after phenotyping

  # 178

  # Gene scoring
  #

}
