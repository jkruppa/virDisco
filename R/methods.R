setMethod(
  f = "[",
  signature = "program_list",
  definition = function(x, i, j, drop) {
    if(i == "bowtie2") {return(x@bowtie2)} else {}
    if(i == "bowtie2_build") {return(x@bowtie2_build)} else {}
    if(i == "pauda") {return(x@pauda)} else {}
    if(i == "pauda_build") {return(x@pauda_build)} else {}
    if(i == "samtools") {return(x@samtools)} else {}
    if(i == "trimmomatic") {return(x@trimmomatic)} else {}
    if(i == "star") {return(x@star)} else {}
    if(i == "pandaseq") {return(x@pandaseq)} else {}
    if(i == "ete3") {return(x@ete3)} else {}
    if(i == "blastn") {return(x@blastn)} else {}
    if(i == "seqtk") {return(x@seqtk)} else {}
  }
)


setReplaceMethod(
  f = "[",
  signature = "program_list",
  definition = function(x, i, j , value) {
    if(i == "bowtie2") {x@bowtie2 <- value} else {}
    if(i == "pauda") {x@pauda <- value} else {}
    if(i == "star") {x@star <- value} else {}
    if(i == "bowtie2_build") {x@bowtie2_build <- value} else {}
    if(i == "pauda_build") {x@pauda_build <- value} else {}
    if(i == "samtools") {x@samtools <- value} else {}
    if(i == "trimmomatic") {x@trimmomatic <- value} else {}
    if(i == "pandaseq") {x@pandaseq <- value} else {}
    if(i == "ete3") {x@ete3 <- value} else {}
    if(i == "blastn") {x@blastn <- value} else {}
    if(i == "seqtk") {x@seqtk <- value} else {}
  }
)


setMethod(
  f = "show",
  signature = "program_list",
  definition = function(object) {
    cat(str_c("Bowtie2 path set to ", object@bowtie2, "\n"))   
    cat(str_c("Pauda path set to ", object@pauda, "\n"))
    cat(str_c("Star path set to ", object@star, "\n"))   
    cat(str_c("Bowtie2_build path set to ", object@bowtie2_build, "\n"))   
    cat(str_c("Pauda_build path set to ", object@pauda_build, "\n"))
    cat(str_c("Samtools path set to ", object@samtools, "\n"))    
    cat(str_c("Trimmomatic path set to ", object@trimmomatic, "\n"))
    cat(str_c("Pandaseq path set to ", object@pandaseq, "\n"))
    cat(str_c("Ete3 path set to ", object@ete3, "\n"))    
    cat(str_c("Blastn path set to ", object@blastn, "\n"))
    cat(str_c("Seqtk path set to ", object@seqtk, "\n"))
  }
)


setMethod(
  f = "[",
  signature = "par_list",
  definition = function(x, i, j, drop) {
    if(i == "genome_index") {return(x@index_genome_dir)} else {}
    if(i == "amino_index") {return(x@index_amino_dir)} else {}
    if(i == "paired") {return(x@paired)} else {}
    if(i == "tax") {return(x@tax)} else {}
    if(i == "qc") {return(x@qc)} else {}
    if(i == "plot") {return(x@plot)} else {}
    if(i == "num_plot") {return(x@num_plot)} else {}
    if(i == "black_list") {return(x@black_list)} else {}
    if(i == "run_dna_mapping") {return(x@run_dna_mapping)} else {}
    if(i == "run_pep_mapping") {return(x@run_pep_mapping)} else {}
    if(i == "consensus") {return(x@consensus)} else {}
    if(i == "gen_prot_sep") {return(x@gen_prot_sep)} else {}
    if(i == "check_host") {return(x@check_host)} else {}
    if(i == "clean") {return(x@clean)} else {}
    if(i == "tmp_dir") {return(x@tmp_dir)} else {}
    if(i == "mapper") {return(x@mapper)} else {}
    if(i == "plot_id") {return(x@plot_id)} else {}
    if(i == "sql_dir") {return(x@sql_dir)} else {}
    if(i == "pdf_file") {return(x@pdf_file)} else {}
    if(i == "ref_seq_file") {return(x@ref_seq_file)} else {}
    if(i == "ncore") {return(x@ncore)} else {}
    if(i == "species_info") {return(x@species_info)} else {}
    if(i == "prot_info") {return(x@prot_info)} else {}
    if(i == "map_dna_in") {return(x@map_dna_in)} else {}
    if(i == "map_pep_in") {return(x@map_pep_in)} else {}
    if(i == "decoy") {return(x@decoy)} else {}
 }
)

setReplaceMethod(
  f = "[",
  signature = "par_list",
  definition = function(x, i, j , value) {
    if(i == "paired") {x@paired <- value} else {}
    if(i == "tax") {x@tax <- value} else {}
    if(i == "qc") {x@qc <- value} else {}
    if(i == "plot") {x@plot <- value} else {}
    if(i == "num_plot") {x@num_plot <- value} else {}
    if(i == "black_list") {x@black_list <- value} else {}
    if(i == "plot_id") {x@plot_id <- value} else {}
    if(i == "run_dna_mapping") {x@run_dna_mapping <- value} else {}
    if(i == "run_pep_mapping") {x@run_pep_mapping <- value} else {}
    if(i == "consensus") {x@consensus <- value} else {}
    if(i == "gen_prot_sep") {x@gen_prot_sep <- value} else {}
    if(i == "check_host") {x@check_host <- value} else {}    
    if(i == "clean") {x@clean <- value} else {}
    if(i == "tmp_dir") {x@tmp_dir <- value} else {}
    if(i == "mapper") {x@mapper <- value} else {}
    if(i == "ncore") {x@ncore <- value} else {}
    if(i == "pdf_file") {x@pdf_file <- value} else {}
    if(i == "ref_seq_file") {x@ref_seq_file <- value} else {}
    if(i == "map_dna_in") {x@map_dna_in <- value} else {}
    if(i == "map_pep_in") {x@map_pep_in <- value} else {}
    if(i == "decoy") {x@decoy <- value} else {}
    return(x)
  }
)

setMethod(
  f = "show",
  signature = "par_list",
  definition = function(object) {
    cat("*** Class par_list, method Show *** \n")
    cat(str_c("mapper = ", object@mapper, "\n"))   
    cat("** Flags are set to: \n")
    cat(str_c("qc = ", object@qc, "\n"))
    cat(str_c("consensus = ", object@consensus, "\n"))    
    cat(str_c("clean = ", object@clean, "\n"))
    cat(str_c("tax = ", object@tax, "\n"))
    cat(str_c("gen_prot_sep = ", object@gen_prot_sep, "\n")) 
    cat("** File paths are set to: \n")
    cat(str_c("pdf_file = '", object@pdf_file, "'\n"))
    cat(str_c("ref_seq_file = '", object@ref_seq_file, "'\n"))
    cat(str_c("tmp_dir = '", object@tmp_dir, "'\n"))
    cat(str_c("sql_dir = '", object@sql_dir, "'\n"))
    cat("** Indexes are set to: \n")
    cat(str_c("index_genome_dir = '", object@index_genome_dir, "'\n"))
    cat(str_c("index_amino_dir = '", object@index_amino_dir, "'\n"))
    cat("** Gene information are: \n")
    if(any(class(par_list@prot_info) == "empty")) {
      cat("No gene information are given\n")
    } else {
      print(head(object@prot_info, n = 1))
    }
    cat("** Species information are: \n")
    if(any(class(par_list@species_info) == "empty")) {
      cat("No species information are given\n")
    } else {
      print(head(object@species_info, n = 1))
    }
  }
)
##' Desc
##'
##' Detail
##' @title Test 
##' @param object 
##' @return boolean
##' @author Jochen Kruppa
##' @export
setGeneric("tbl_sqlite_check", function(object) {standardGeneric("tbl_sqlite_check")})

setMethod(
  f = "tbl_sqlite_check",
  signature = "par_list",
  definition = function(object) {
    check_aa_info <- FALSE
    check_species_info <- FALSE
    ## check the names of the info files
    names_aa_info <- c("prot_id", "genebank_id", "pos_start", "pos_end")
    if(!all(names_aa_info %in% suppressWarnings(names(collect(object@prot_info, n = 1))))){
      stop("prot_info must include the following columns: ", str_c(names_aa_info, collapse = ", "), ". Use setup_aa_info_sqlite() to setup the right sqlite database")
    } else {
      check_aa_info <- TRUE
    }
    names_species_info <- c("genebank_id", "Length", "Description")
    if(!all(names_species_info %in% suppressWarnings(names(collect(object@species_info, n = 1))))){
      stop("prot_info must include the following columns: ", str_c(species_info, collapse = ", "), ". Use setup_species_info_sqlite() to setup the right sqlite database")
    } else {
      check_species_info <- TRUE
    }
    if(all(check_species_info, check_aa_info)){
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
)


    

