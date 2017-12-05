
setClass(
  Class = "program_list",
  representation = representation(
    bowtie2 = "character",
    bowtie2_build = "character",
    star = "character",
    pauda = "character",
    pauda_build = "character",    
    pandaseq = "character",
    ete3 = "character",
    blastn = "character",
    seqtk = "character",
    samtools = "character",
    trimmomatic = "character"
  ),
  validity = function(object){
    if(!file.exists(object@ete3)) {
      stop("Cannot find dir 'ete3': No such file or directory")
    }
    if(!file.exists(object@blastn)) {
      stop("Cannot find dir 'blastn': No such file or directory")
    }
    if(!file.exists(object@seqtk)) {
      stop("Cannot find dir 'seqtk': No such file or directory")
    }
    if(!file.exists(object@pandaseq)) {
      stop("Cannot find dir 'pandaseq': No such file or directory")
    }
    if(!file.exists(object@star)) {
      stop("Cannot find dir 'star': No such file or directory")
    }
    if(!file.exists(object@bowtie2)) {
      stop("Cannot find dir 'bowtie2': No such file or directory")
    }
    if(!file.exists(object@pauda)) {
      stop("Cannot find dir 'pauda': No such file or directory")
    }
    if(!file.exists(object@pauda_build)) {
      stop("Cannot find dir 'pauda': No such file or directory")
    }
    if(!file.exists(object@bowtie2_build)) {
      stop("Cannot find dir 'pauda': No such file or directory")
    }
    if(!file.exists(object@samtools)) {
      stop("Cannot find dir 'samtools': No such file or directory")
    }
    if(!file.exists(object@trimmomatic)) {
      stop("Cannot find dir 'trimmomatic': No such file or directory")
    }
    return(TRUE)
  }
)


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
    ncore = "numeric",
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
    ncore = 1,
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

