#' Perform a descendant search for a query term
#'
#' This function will take a query term, and find all descendants and neighbors of that term
#' in the specified ontologies. By default, a match is considered if a the query is a subset
#' of a database entry. i.e., "alzheimer" matches "Alzheimer's disease".
#' This process is important for increasing the chances of
#' getting a string match in the disease-gene databases.
#'
#' @param term A query term
#' @param ontologies Defaults to `c("ctd", "doid")`
#' @param exact_match Defaults to `FALSE`. If true, no partial string mathching is performed.
#'
#' @seealso \code{\link{output_gene_prioritization}}, where this function is called
#'
#'
#' @export
disease_extension <- function(term, ontologies = c("ctd", "doid"), exact_match = FALSE) {

  # If the ontology files aren't in the global environment, get them
  if(!(exists("DOID") & exists("CTD"))) retrieve_ontologies()

  # Idiot check
  if(!("ctd" %in% ontologies | "doid" %in% ontologies)) stop("Please specify either ctd or doid")

  term <- term

  if("doid" %in% ontologies){
    doid_results <- .ontology_search(term, DOID)

    message("Finished searching the Disease Ontology Database")
    if (length(doid_results > 0)) {
      doid_results <-
        cbind(doid_results, rep("DISEASE_ONTOLOGY", length(doid_results)))
    }
  }
    else
      doid_results <-
      tibble(a = character(0), b = character(0)) # If no results found, make a dummy variable


  if("ctd" %in% ontologies){
    ctd_results <- .ontology_search(term, CTD)

    message("Finished searching the CTD Medic Database")
    if (length(ctd_results) > 0) {
      ctd_results <-
        cbind(ctd_results, rep("CTD_DISEASE", length(ctd_results)))
    }
  }
  else
    ctd_results <-
    tibble(a = character(0), b = character(0))  # If no results found, make a dummy variable

  output <- as_tibble(rbind(doid_results, ctd_results))
  if(length(output) == 0) stop("No matching terms were found for you query :-( Please check your spelling!")

  # Include the original input in the output table. This can help when the actual term isn't in the database
  # i.e. "lymphoma" has a lot of descendants, is no in the DB by itself
  output <- rbind(tibble(doid_results = term, V2 = "MY_TERM"), output)
  names(output) <- c("Disease", "Source")
  output$Source <- output$Source %>% as.factor
  output

}

#' Retrieve genes related to a disease query
#'
#' This function calls `disease_extension` to perform a descendant search, then prioritizes the
#' disease terms matched to the disease-gene databases using system described in the paper.
#'
#' In brief, each database has its own quality metrics. These are given an ad hoc score, and matches
#' are weighted based on this score. The prioritized list is ordered, and normalized by diving through
#' all scores by the max score.
#' @param disease_term A disease for which to find related genes
#' @param ontologies Defaults to c("ctd", "doid").
#' @param return_output Logical. Specify if the output should be returned. Useful if you only want a wordcloud.
#' @param output_mode Specifies if the output is a table of genes, or a list of tibbles with individual genes and source information
#'
#' @examples
#' output_gene_prioritization("alzheimer's disease")
#' output_gene_prioritization("alzheimer's disease", ontologies = "doid")

#' @export
output_gene_prioritization <- function(disease_term,
                                       ontologies = c("ctd", "doid"),
                                       output_mode = "aggregate",
                                       return_output = TRUE,
                                       wordcloud = FALSE,
                                       predict = FALSE)
  {

  useless_words <- c("disease", "syndrome")

  term_k <- disease_term %>% tolower %>% .cleaner



  # Make sure databases are all there
  if(!(exists("DOID") & exists("CTD"))) build_gene_id_syn()

  # Run extension on each Term_k
  diseases_k <- lapply(term_k, function(k) disease_extension(k, ontologies = ontologies))
  if (max(sapply(diseases_k, nrow)) == 0)
    stop("No matching disease terms found :-( Please check your spelling!!")

  message("Disease term extension completed")

  # Clean each Term_k and remove duplicates
  for (i in 1:length(diseases_k)) {
    diseases_k[[i]]$Disease %<>% .cleaner
  }
  diseases_k <-
    lapply(diseases_k, function(x)
      x %>% mutate(Disease = Disease %>% tolower %>% trimws) %>% distinct)
  names(diseases_k) <- term_k

  # Make sure DB_COMPILED_DISEASE_SCORES is in the global environnment
  if(!exists("DB_COMPILED_GENE_DISEASE_SCORES")) build_gene_disease_reference()

  # Match each disease in diseases_k with the precompiled database
  #### the most important part of this whole thing!!! ####

  message("Matching to database with extended terms!")
  disease_scores1 <- lapply(diseases_k, function(x) {
    DB_COMPILED_GENE_DISEASE_SCORES %>%
      filter( tolower(.cleaner(.$DISEASE)) %>% map_lgl( ~ . %in% x$Disease) )
    })

  # Match to the precompiled database by greping the original term to the precompiled database

  message("Matching to database with input term!")
  disease_scores2 <- lapply(term_k, function(k) {
    DB_COMPILED_GENE_DISEASE_SCORES %>%
      filter(grepl(k, tolower(.cleaner(DISEASE)), ignore.case = TRUE))
  })


  # Combine these results and take the distinct matches
  disease_scores <- c(disease_scores1, disease_scores2)

  if (max(sapply(disease_scores, nrow)) == 0)
    stop("No disease / gene relationships found. Please adjust your query!!")

  # Aggregate the disease_k list

  if (output_mode == "aggregate"){
    out <-
      # For each disease in diseases_k:
      # aggregate all genes and sum scores of genes which appear more than once
      # divide through by the total genes that appear
      lapply(disease_scores, function(disease_k)
      {disease_k %>% group_by(GENE) %>% summarise(total_SCORE = sum(SCORE) / nrow(disease_k))}) %>%

      # Combine all k tibbles
      # This is the first time the Term_k's from the input are combined!!!
      do.call("rbind", .) %>%

      # Aggregate by gene again
      group_by(GENE) %>% summarise(reported_SCORE = sum(total_SCORE)) %>%

      # Sort the big tibble by decreasing score
      arrange(desc(reported_SCORE))

    # Divide scores by the max (normalization step)
    out$reported_SCORE <- out$reported_SCORE / out$reported_SCORE[1]

    # Add gene numbers to the final output tibble
    out <- left_join(out, HUMAN_GENE_ID %>% select(GeneID, Symbol), by = c("GENE" = "Symbol"))

    # Get it in the right order
    out %<>% select("GeneID", "GENE", "reported_SCORE")
  } else if (output_mode == "list"){
    out <-
      disease_scores %>% do.call("rbind", .) %>% split(.$GENE)

    out_score <- sapply(out, function(x) x$SCORE %>% sum) %>% sort(decreasing = TRUE)
    out_score_norm <- out_score / max(out_score) %>% round(3)

    out <- out[names(out_score)]
    names(out) <- paste(names(out_score), "Raw score:", out_score, "Normalized score:", out_score_norm)
  }

  if(wordcloud){
    pal <- brewer.pal(6, "Dark2")
    par(mar = rep(0,4))
    wordcloud(
      out$GENE,
      out$reported_SCORE * 500,
      colors = pal,
      max.words = 100
    )
  }

  if(predict){
    message("The gene prediction process starts")
    output_predicted_genes(out)
  } else if(return_output) return(out)

}

#' Write an itemized gene list to file
#'
#' If you want to see the sources for each gene call, this function will write the results
#' of `output_gene_prioritization` with `output_mode == "list"` to a text file that is conveniently
#' readable.
#'
#' @param itemized_gene_list The output of list `output_gene_prioritization` with `output_mode == "list"`
#' @param filename File name to write to the the filesystem
#'
#' @export
write_gene_list <- function(itemized_gene_list, filename){

  # Write the file header
  write_lines(paste("Number of genes identified in the databases:", length(itemized_gene_list)), filename)

  # For each gene in the gene list, write the gene rank information, then write the table
  for (gene in names(itemized_gene_list)) {
    write_lines(gene, filename, append = TRUE)
    write_tsv(as.data.frame(itemized_gene_list[gene]), filename, append = TRUE)
    write_lines("", filename, append = TRUE)
  }
}
