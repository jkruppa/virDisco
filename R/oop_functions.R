##' Test
##'
##' Test
##' @title TEst 
##' @param index_genome_dir 
##' @param index_amino_dir 
##' @param qc 
##' @param set_mapper 
##' @param ref_seq_file 
##' @param prot_info_df 
##' @param species_info_df 
##' @return foo
##' @author Jochen Kruppa
##' @export
set_par_list <- function(index_genome_dir,
                         index_amino_dir,
                         qc = TRUE,
                         set_mapper = "bowtie",
                         ref_seq_file = "",
                         prot_info_df = "",
                         species_info_df = "")
{
  if(any(nchar(prot_info_df) == 0) | any(nchar(species_info_df) == 0)) {
    new(Class = "par_list",
        index_genome_dir = index_genome_dir,
        index_amino_dir = index_amino_dir,
        ref_seq_file = ref_seq_file,
        mapper = set_mapper,
        qc = qc)
  } else {
    new(Class = "par_list",
        index_genome_dir = index_genome_dir,
        index_amino_dir = index_amino_dir,
        qc = qc,
        mapper = set_mapper,
        ref_seq_file = ref_seq_file,
        prot_info = prot_info_df,
        species_info = species_info_df)
  }
}

##' Test
##'
##' Test
##' @title Test 
##' @param bowtie_dir 
##' @param pauda_dir 
##' @param samtools_dir 
##' @param trimmomatic_dir 
##' @return program_list
##' @author Jochen Kruppa
##' @export
set_program_list <- function(bowtie_dir,
                             pauda_dir,
                             samtools_dir,
                             trimmomatic_dir)
{
  new(Class = "program_list",
      bowtie = bowtie_dir,
      pauda = pauda_dir,
      samtools = samtools_dir,
      trimmomatic = trimmomatic_dir)
}
