##' The function generates the SQlite databse needed for the plotting.
##'
##' The function generates the SQlite databse needed for the plotting
##' of the hit seqeunces.
##' @title Function to build up the aa_info_sqlite database
##' @param ncbi_aa_seq A AAStringSet containing amino acids
##' @param db_file File path to the sqlite3 database
##' @return Null
##' @author Jochen Kruppa
##' @export
setup_aa_info_sqlite <- function(ncbi_aa_seq, db_file){
  talk("Build info data frame")
  if(file.exists(db_file)) stop(str_c(db_file, " exists and must be deleted first!"))
  aa_pos_df <- ldply(names(ncbi_aa_seq), function(x){
    if(nchar(x) > 40) {
      info <- str_split(str_split(x, " ", simplify = TRUE)[,1], "_", simplify = TRUE)
      return(tibble(prot_id = info[1], genebank_id = info[1],
                    pos_start = info[2], pos_end = info[3]))
    } else {
      tmp_df <- tbl_df(str_split(x, "_", simplify = TRUE))
      if(length(tmp_df) == 5){
        tmp_df <- tibble(tmp_df[[1]], str_c(tmp_df[2], "_", tmp_df[3]), tmp_df[[4]], tmp_df[[5]])
      }
      names(tmp_df) <- c("prot_id", "genebank_id", "pos_start", "pos_end")
      return(tmp_df)
    }
  }, .parallel = FALSE, .progress = "text")
  aa_pos_df$ind <- 1:nrow(aa_pos_df)
  aa_pos_df$pos_start <- as.numeric(aa_pos_df$pos_start)  
  aa_pos_df$pos_end <- as.numeric(aa_pos_df$pos_end)
  ## build up sql lite data base
  talk("Setup SQLite")
  aa_info_db <- src_sqlite(db_file, create = TRUE)
  aa_info_sqlite <- copy_to(aa_info_db, aa_pos_df, temporary = FALSE)
}



##' The function generates the SQlite databse needed for the plotting.
##'
##' The function generates the SQlite databse needed for the plotting
##' of the hit seqeunces.
##' @title Function to build up the species_info_sqlite database
##' @param genebank_id Vector of genebank_ids
##' @param length Length vector of the sequence
##' @param description Description vector of the genebank_ids
##' @param tax_id NCBI tax_id connected with the genebank_ids
##' @param mol_type Molecular type
##' @param strain Strain
##' @param accession Accession number
##' @param organism Organsim
##' @param db_file File path to save the sqlite3 database
##' @return NULL
##' @author Jochen Kruppa
##' @export
setup_species_info_sqlite <- function(genebank_id,
                                      length,
                                      description,
                                      tax_id = NA,
                                      mol_type = NA,
                                      strain = NA,
                                      accession = NA,
                                      organism = NA,
                                      db_file){
  talk("Setup db at ", db_file)
  if(file.exists(db_file)) stop(str_c(db_file, " exists and must be deleted first!"))
  species_info_df <- tibble(genebank_id,
                            mol_type,
                            strain,
                            tax_id,
                            Accession = accession,
                            Length = length,
                            Organism = organism,
                            Description = description)
  ## build up sql lite data base
  species_info_db <- src_sqlite(db_file, create = TRUE)
  species_info <- copy_to(species_info_db, species_info_df, temporary = FALSE)
}

##' The function generates the SQlite databse needed for the plotting.
##'
##' The function generates the SQlite databse needed for the plotting
##' of the hit seqeunces.
##' @title Function to build up the aa_info_sample_sqlite database
##' @param prot_id Protein ID of the amino acids
##' @param genebank_id Genebank ID connected to the protein IDs
##' @param prot_start Protein start position on the sequence
##' @param prot_end Protein end position on the sequence
##' @param db_file File path to save the sqlite3 database
##' @return NULL
##' @author Jochen Kruppa
##' @export
setup_aa_info_sample_sqlite <- function(prot_id, genebank_id, prot_start, prot_end, db_file){
  talk("Setup db at ", db_file)
  if(file.exists(db_file)) stop(str_c(db_file, " exists and must be deleted first!"))
  aa_info_df <- tibble(prot_id,
                       genebank_id,
                       pos_start = prot_start,
                       pos_end = prot_end)
  ## build up sql lite data base
  aa_info_db <- src_sqlite(db_file, create = TRUE)
  aa_info <- copy_to(aa_info_db, aa_info_df, temporary = FALSE)
}
