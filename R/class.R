## we need to define the class tbl_sqlite
setOldClass("tbl_dbi")

## define class par_list
setClass(
  Class = "par_list",
  representation = representation(
    tax = "logical",
    qc = "logical",
    clean = "logical",
    paired = "logical",
    plot = "logical",
    black_list = "character",
    gen_prot_sep = "character",
    check_host = "logical",
    run_dna_mapping = "logical",
    run_pep_mapping = "logical",
    consensus = "logical",
    tmp_dir = "character",
    sql_dir = "character",
    index_genome_dir = "character",
    index_amino_dir = "character",
    ref_seq_file = "character",
    prot_info = "tbl_dbi",
    species_info = "tbl_dbi",
    plot_id = "character",
    mapper = "character",
    num_plot = "numeric",
    pdf_file = "character"
  ),
  prototype = prototype(
    tax = TRUE,
    clean = TRUE,
    paired = FALSE,
    plot = TRUE,
    plot_id = "",
    run_dna_mapping = TRUE,
    run_pep_mapping = TRUE,
    consensus = TRUE,
    gen_prot_sep = "_",
    black_list = "",
    check_host = TRUE,
    mapper = "bowtie",
    num_plot = 25,
    prot_info = structure(tibble(), class = c("tbl_dbi", "empty")),
    species_info = structure(tibble(), class = c("tbl_dbi", "empty")),
    tmp_dir = file.path("/home/temp"),
    sql_dir = file.path("/home/sql")
  ),
  validity = function(object){
    if(nchar(object@index_genome_dir) == 0) {
      stop("index_genome_dir must be given!")
    }
    if(object@ref_seq_file != ""){
      if(!file.exists(object@ref_seq_file)) {
        stop("cannot find dir 'ref_seq_file': No such file or directory")
      }
    }
    if(!file.exists(dirname(object@index_genome_dir))) {
      stop("cannot find dir 'index_genome_dir': No such file or directory")
    }
    if(!file.exists(object@index_amino_dir)) {
      stop("cannot find dir 'index_amino_dir': No such file or directory")
    }
    if(!file.exists(object@tmp_dir)){
      stop("cannot find dir 'tmp_dir': No such file or directory")
    }
    if(!file.exists(object@sql_dir)) {
      stop("cannot find dir 'sql_dir': No such file or directory")
    }
    return(TRUE)
  }
)

