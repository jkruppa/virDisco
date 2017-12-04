##' Test
##'
##' Test
##' @title Test 
##' @param dna_set 
##' @param aa_set 
##' @param index_dir 
##' @param index_name 
##' @param bowtie 
##' @param pauda 
##' @param star 
##' @param force 
##' @param tmp_dir 
##' @return file.list 
##' @author Jochen Kruppa
##' @export
build_index <- function(dna_set, aa_set,
                        index_dir = NULL,
                        index_name = NULL,
                        bowtie = TRUE,
                        pauda = TRUE,
                        star = FALSE,
                        force = FALSE,
                        tmp_dir = tempdir()){
  talk("Setup files")
  ## file handling for the analysis
  if(is.null(index_dir)) {
    stop("No index_dir is given!")
  }
  if(is.null(index_name)) {
    stop("No index_name is given!")
  }
  if(!file.exists(index_dir)) dir.create(index_dir)
  ## check for file or stringset
  if(class(dna_set) == "DNAStringSet") {
    dna_set_fa <- file.path(index_dir, "dna_set.fa")
    writeXStringSet(dna_set, dna_set_fa)
  } else {
    dna_set_fa <- dna_set
  }
  if(class(aa_set) == "AAStringSet") {
    aa_set_fa <- file.path(index_dir, "aa_set.fa")
    writeXStringSet(aa_set, aa_set_fa)
  } else {
    aa_set_fa <- aa_set
  }
  ## generate index
  star_index <- file.path(index_dir, "star")
  bowtie_index <- file.path(index_dir, "bowtie", str_c(index_name, "_bowtie"))
  pauda_tmp_index <- file.path(tmp_dir, "pauda")
  ## file checking
  if(file.exists(dirname(bowtie_index)) & !force & bowtie){
    bowtie <- FALSE
    warning("Bowtie index exists. Skiping. Use 'force = TRUE' to overwrite...")
  }
  if(file.exists(dirname(star_index)) & !force & star){
    star <- FALSE
    warning("Star index exists. Skiping. Use 'force = TRUE' to overwrite...")
  }
  if(file.exists(dirname(pauda_tmp_index)) & !force & pauda){
    pauda <- FALSE
    warning("Pauda index exists. Skiping. Use 'force = TRUE' to overwrite...")
  }
  ## start bowtie
  if(bowtie){
    if(!file.exists(dirname(bowtie_index))) dir.create(dirname(bowtie_index))
    talk("Start bowtie build")
    bowtie2BuildCMD <- paste(program_list["bowtie2_build"],
                             "-q",
                             dna_set_fa,
                             bowtie_index,
                             "1>", file.path(tmpDir, "bowtie.log"), "2>&1")
    runCMD(bowtie2BuildCMD)
  }
  ## start star
  if(star){
    talk("Start star build")
    if(!file.exists(star_index)) dir.create(star_index)
    star_build_CMD <- paste(STAR,
                            "--runMode", "genomeGenerate",
                            "--genomeDir", star_index, 
                            "--genomeFastaFiles", dna_set_fa,
                            "--runThreadN", par_list["ncore"],
                            "1>", file.path(tmp_dir, "mapper.log"), "2>&1")
    runCMD(star_build_CMD)
  }  
  if(pauda){ 
    ## build PAUDA index by bowtie2
    if(file.exists(pauda_tmp_index)) unlink(pauda_tmp_index, recursive = TRUE)
    ## do NOT create the dir...
    talk("Start pauda build")
    runCMD(paste(program_list["pauda_build"],
                 aa_set_fa,
                 pauda_tmp_index,
                 "1>", file.path(tmp_dir, "pauda.log"), "2>&1"))
    pauda_files <- list.files(tmp_dir, pattern = "pauda", full.names = TRUE,
                              include.dirs = TRUE)
    file.copy(pauda_files, index_dir, recursive = TRUE)
    ## and remove them from working dir...
    unlink(pauda_files, recursive = TRUE)
  }
  return(list(bowtie_index = bowtie_index,
              star_index = star_index,
              pauda_index = file.path(index_dir, "pauda")))
}
