### Right now, I'm just going to focus on extending the disease using terms in disease ontology
# and CTD Disease. I think this should cover the bases pretty well.


#' Retrieve the obo ontology files needed for disease term descendant searching.
#' Currently, this package loads the CTD Disease ontology and the Disease ongology (doid.obo)
#'
retrieve_ontologies <- function(){
  message("Downloading the Disease ontology")
  doid <-
    get_ontology("http://purl.obolibrary.org/obo/doid.obo", extract_tags = "everything")

  message("Downloading the Clinical Toxicogenomics Database - Disease Names")
  temp <- tempfile()
  download.file("http://ctdbase.org/reports/CTD_diseases.obo.gz", temp)
  ctd <-
    get_ontology(temp, extract_tags = "everything") #Pause on this one for now...
#  hpo <-
#    get_ontology("http://purl.obolibrary.org/obo/hp.obo", extract_tags = "everything")
}

cleaner <- function(query) {
  query %>%
    str_remove("'s") %>%
    str_remove("[Ll]ate onset|[Ff]amilial") %>%
    str_remove(" \\d+$")
}

ontology_search <- function(term, ontology, exact_match = FALSE) {
  if (class(ontology) != "ontology_index" |
      missing(ontology))
    stop("Please supply a valid ontology file")

  term %<>% cleaner
  db_names_clean <- cleaner(ontology$name)
  db_syn_clean <- cleaner(ontology$synonym)

  if (exact_match) {
    id <-
      ontology$id[grep(
        db_names_clean,
        pattern = paste("^", term, "$", sep = ""),
        ignore.case = TRUE
      )]
  }
  else{
    id <-
      ontology$id[grep(db_names_clean,
                       pattern = term,
                       ignore.case = TRUE)]
    syn <-
      ontology$id[grep(db_syn_clean, pattern = term, ignore.case = TRUE)]
    id <- union(id, syn)
  }

  # Expand the list by adding children iteratively unil the list doesn't grow any more
  id_length <- length(id) - 1
  while (length(id) > id_length) {
    id_length <- length(id)
    id <- c(id, ontology$children[id]) %>% unlist %>% unique
  }

  # Get names of all terms and ddd synonym terms to the names list
  names <- ontology$name[id]
  names <-
    c(names,
      ontology$synonym[id] %>% unlist %>% gsub("^\\\"", "", .) %>% gsub("\\\".+", "", .)) # add synonym terms
  names
}



disease_extension <- function(term, exact_match = FALSE) {
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

