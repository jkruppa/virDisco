##' Test
##'
##' Test
##' @title Test 
##' @param genebank_id 
##' @param length 
##' @param description 
##' @param tax_id 
##' @param mol_type 
##' @param strain 
##' @param accession 
##' @param organism 
##' @param db_file 
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

##' Test
##'
##' Test
##' @title Test 
##' @param prot_id 
##' @param genebank_id 
##' @param prot_start 
##' @param prot_end 
##' @param db_file 
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
