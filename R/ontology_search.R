#' Retrieve the obo ontology files needed for disease term descendant searching.
#' Currently, this package loads the CTD Disease ontology and the Disease ongology (doid.obo)
#' @export

retrieve_ontologies <- function(){
  message("Downloading the Disease ontology")
  doid <<-
    get_ontology("http://purl.obolibrary.org/obo/doid.obo", extract_tags = "everything")

  message("Downloading the Clinical Toxicogenomics Database - Disease Names")
  temp <- tempfile()
  download.file("http://ctdbase.org/reports/CTD_diseases.obo.gz", temp)
  ctd <<-
    get_ontology(temp, extract_tags = "everything") #Pause on this one for now...
#  hpo <-
#    get_ontology("http://purl.obolibrary.org/obo/hp.obo", extract_tags = "everything")
}

#'
.cleaner <- function(query) {
  query %>%
    str_remove("'s") %>%
    str_remove("[Ll]ate onset|[Ff]amilial") %>%
    str_remove(" \\d+$")
}

#' @export
ontology_search <- function(term, ontology, exact_match = FALSE) {
  if (class(ontology) != "ontology_index" |
      missing(ontology))
    stop("Please supply a valid ontology file")

  term %<>% .cleaner
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
