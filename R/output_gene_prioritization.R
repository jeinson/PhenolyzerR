useless_words <- c("disease", "syndrome")

output_gene_prioritization <- function(disease_term) {

  # Split the query string into individual words and get rid of common "useless words"
  term_k <-
    disease_term %>% str_split(" ") %>% unlist %>% discard( ~ .x %in% useless_words)

  if (length(term_k) > 1){print(paste("Splitting your query into", length(term_k), "searches"))}

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

  # Match each disease in diseases_k with the precompiled database
  disease_scores <- lapply(diseases_k, function(x) {
    DB_COMPILED_GENE_DISEASE_SCORES %>% filter(tolower(.$DISEASE) %>% map_lgl( ~ . %in% x$Disease))
  })

  if (max(sapply(relevant_disease_scores, nrow)) == 0)
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