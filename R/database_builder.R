#' This function downloads a list of all genes and saves a tibble of Entrez gene IDs, gene number,
#' and synonyms in the global environment
#' @example build_gene_id_syn() -> GENE_ID
#' @export
build_gene_id_syn <- function(){
  #### RefSeq Gene IDs ####
  message("Downloading the RefSeq gene information table")
  suppressMessages(
    HUMAN_GENE_ID <<- read_tsv("ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz", progress = TRUE) %>%
      select(GeneID, Symbol, Synonyms)
  )
}

#' Aggregate all the databases needed to run PhenolyzerR
#'
#' `build_gene_disease_reference` is downloads gene-disease mappings from OMIM, ClinVar, and the GWAS repository.
#' This will only need to be run once. When the script finishes running,
#' you will have a table called `DB_COMPILED_GENE_DISEASE_SCORES` in your local environment.
#'
#' @return A Tibble
#'
#' @export
build_gene_disease_reference <- function(){
  # Compile the Gene - Disease Scores
  # This is the meat of this whole project!!

  #### OMIM ####
  # Get OMIM Genemap
  # Hope I'm not breaking any laws by doing this...
  message("Downloading the OMIM genemap")
  tryCatch(
    suppressMessages(suppressWarnings(
      OMIM_conf <- read_tsv("https://data.omim.org/downloads/WIFQrs35Rw-0RbsyXS9pWA/genemap.txt", skip = 3) %>%
        select("Gene Symbols", "MIM Number", "Confidence")
    )),
    error = function(e){message("Unable to Connect to OMIM. Please try again!")}
  )
    suppressMessages(suppressWarnings(
    MIMs <- read_tsv("https://data.omim.org/downloads/WIFQrs35Rw-0RbsyXS9pWA/morbidmap.txt", skip = 3) %>% # Download the morbidmap
      select(-`Cyto Location`, -`Gene Symbols`) # Get rid of cyto location and gene symbols column
    ))
    OMIM <- left_join(OMIM_conf, MIMs, "MIM Number") # Join these tables together



  S_SCORE_OMIM <- list("C" = 1, "P" = 0.75, "L" = 0.5, "I" = 0.25)

  # These function are used in the next step::
  cleanup <- function(term){
    term %>%
      str_remove("\\(\\d\\)") %>%
      str_remove(", \\d+ $") %>%
      gsub("^\\W(.*?)\\W*susc?eptibi?lity( to)?,?.*", "\\1", ., perl = TRUE) %>%
      str_replace("[Aa]utism \\d+", "autism") %>%
      str_replace("\\berthematosus\\b", "erythematosus")
  }

  eliminateNonwords <- function(term){
    term %>%
      str_remove_all(regex("[[:punct:]]")) %>%
      str_remove("\\bwith\\b.*$") %>%
      str_remove("\\bwithout\\b.*$") %>%
      trimws
  }

  OMIM_SCORES <- OMIM %>% select("Gene Symbols", "# Phenotype", "MIM Number", "Confidence")
  OMIM_SCORES <- OMIM_SCORES[complete.cases(OMIM_SCORES),] %>% # Remove NA Rows

    # Make sure all confidence codes are capital letters
    mutate(Confidence = toupper(Confidence)) %>%

    # Get rid of rows which have something other than the 4 allowed codes
    filter(Confidence %in% c("C", "P", "L", "I")) %>%

    # Convert confidence codes to scores
    mutate("Score" = unlist(S_SCORE_OMIM[Confidence])) %>%

    # Add a Source column
    mutate("Source" = "OMIM") %>%

    # Get rid of the confidence codes
    select(-Confidence) %>%

    # Add "OMIM" ot source column
    mutate(`MIM Number` = paste("OMIM:", `MIM Number`, sep = "")) %>%

    # Get rid of comma separation in gene column, so every column has only one gene
    mutate(`Gene Symbols` = strsplit(as.character(`Gene Symbols`), ", ")) %>%
    unnest(`Gene Symbols`) %>%

    # Cleanup the Disease names column
    mutate(`# Phenotype` = `# Phenotype` %>% cleanup %>% eliminateNonwords)

  OMIM_SCORES %<>% select(`Gene Symbols`, `# Phenotype`, `MIM Number`, Score, Source)
  names(OMIM_SCORES) <- c("GENE", "DISEASE", "DISEASE_ID", "SCORE", "SOURCE")


  #### ClinVar ####
  message("Downloading ClinVar gene_condition_source_id file")
  suppressMessages(clinVar <- read_tsv("ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/gene_condition_source_id"))
  clinVar %<>% select(AssociatedGenes, RelatedGenes, DiseaseName, DiseaseMIM)
  clinVar$genes <- paste(clinVar$AssociatedGenes, clinVar$RelatedGenes) %>% gsub("NA", "", .) %>% trimws
  clinVar %<>% select(-AssociatedGenes, -RelatedGenes)
  clinVar <- clinVar[complete.cases(clinVar),]
  clinVar <- aggregate(list(count=rep(1,nrow(clinVar))), clinVar, length) %>% as.tibble

  SCORE_ClinVar <- function(x){
    if (x == 4) return (1.00 * 0.25)
    if (x == 3) return (0.75 * 0.25)
    if (x == 2) return (0.50 * 0.25)
    if (x == 1) return (0.25 * 0.25)
  }

  ClinVar_SCORES <- clinVar %>% select(genes, DiseaseName, DiseaseMIM, count)
  ClinVar_SCORES <- ClinVar_SCORES[complete.cases(ClinVar_SCORES),] %>% # Get rid of missing values
    mutate(Score = sapply(count, SCORE_ClinVar)) %>% # Compute the score of each one based on the alignment
    mutate(DiseaseMIM = paste("OMIM:", DiseaseMIM, sep = "")) %>%
    select(-count) %>% mutate(SOURCE = "ClinVar")
  names(ClinVar_SCORES) <- c("GENE", "DISEASE", "DISEASE_ID", "SCORE", "SOURCE")

  #### GWAS ####
  # Get GWAS associations
  message("Downloading GWAS gene condition associations")
  suppressMessages(suppressWarnings(
    GWAS <- read_tsv("https://www.ebi.ac.uk/gwas/api/search/downloads/full")
    ))

  GWAS_SCORE_WEIGHT <- 0.25

  # This function gets rid of anything in parentheses
  getRidOfAnnotation <- function(x){x %>% str_remove("\\(.*?\\)")}

  GWAS_SCORES <- GWAS %>%
    # Select the columns I want
    select(`REPORTED GENE(S)`, `DISEASE/TRAIT`, PUBMEDID, `P-VALUE`) %>%

    # Get a raw score by calculaitn 1 - p.value and multiple by score weight
    mutate(RAW_SCORE = (1 - `P-VALUE`) * GWAS_SCORE_WEIGHT) %>%
    select(-`P-VALUE`) %>%

    # Remove annotations from disease names
    mutate(`DISEASE/TRAIT` = getRidOfAnnotation(`DISEASE/TRAIT`)) %>%

    # Remove intergenic of unknown genes
    filter(!grepl("[Ii]ntergenic|other|NR", `REPORTED GENE(S)` )) %>%

    # Add "PUBMED" to PUBMEDID
    mutate(PUBMEDID = str_c("PUBMED:", PUBMEDID)) %>%

    # Get rid of duplicate rows and rows with no scores
    distinct %>% filter(!is.na(RAW_SCORE)) %>%

    # Add Column for source
    mutate(SOURCE = "GWAS") %>%

    # Get rid of comma separation in gene column, so every column has only one gene
    mutate(`REPORTED GENE(S)` = strsplit(as.character(`REPORTED GENE(S)`), ", ")) %>%
    unnest(`REPORTED GENE(S)`) %>% select(5,1,2,3,4)

  names(GWAS_SCORES) <- c("GENE", "DISEASE", "DISEASE_ID", "SCORE", "SOURCE")


  #### Put Everything Together ####
  DB_COMPILED_GENE_DISEASE_SCORES <- rbind(OMIM_SCORES, ClinVar_SCORES, GWAS_SCORES)

  #### Fix Synonymous Gene Names from Different DBs ####
  if (!exists("HUMAN_GENE_ID")) build_gene_id_syn()
  syns <- HUMAN_GENE_ID$Synonyms %>% strsplit(split = "\\|")
  names(syns) <- HUMAN_GENE_ID$Symbol
  need_fixin <- DB_COMPILED_GENE_DISEASE_SCORES[DB_COMPILED_GENE_DISEASE_SCORES$GENE %in% unlist(syns),]
  good <- DB_COMPILED_GENE_DISEASE_SCORES[!(DB_COMPILED_GENE_DISEASE_SCORES$GENE %in% unlist(syns)),]

  # Make a Thesaurus of gene syonyms, and fix ones which have a gene synonym as their name
  genes <- sapply(syns, length)
  genes <- rep(names(genes), genes)
  syns %<>% unlist(use.names = FALSE)
  names(genes) <- syns
  need_fixin %<>% mutate(GENE = genes[GENE])

  DB_COMPILED_GENE_DISEASE_SCORES <<- rbind(good, need_fixin) %>% distinct %>% mutate(GENE = trimws(GENE))
}

