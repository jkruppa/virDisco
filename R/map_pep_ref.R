##' Test
##'
##' Test
##' @title Test 
##' @param infile 
##' @param outfile 
##' @param par_list 
##' @return data.frame
##' @author Jochen Kruppa
##' @export
map_pep_ref <- function(infile, outfile, par_list){
  ## mapping to pep reference
  if(par_list["paired"]) {
    infile <- unlist(infile)
    talk("Run pandaseq to build up one fastq file")
    log_file <- file.path(tmpDir, "pandaseq.log")
    panda_cmd <- paste(program_list["pandaseq"],
                       "-f", infile[1],
                       "-r", infile[2],
                       "-F",
                       "-g", log_file,
                       "-w", gsub(".blastx", "_combinded.fastq", outfile))
    runCMD(panda_cmd)
    unlink(log_file)
    ## ------------------------------------------------------------
    ## check for bad index seqeunce    
    fq2fa(gsub(".blastx", "_combinded.fastq", outfile),
          gsub(".blastx", "_combinded.fa", outfile))
    ## read the fasta
    combined_fa <- readDNAStringSet(gsub(".blastx", "_combinded.fa", outfile))
    names(combined_fa) <- str_replace(names(combined_fa), "(^.*):.*", "\\1:1")
    ## write the fasta
    writeXStringSet(combined_fa, gsub(".blastx", "_combinded_mod.fa", outfile))
  } else {
    if(grepl(".gz$", infile)){
      runCMD(paste("gunzip -f", infile))
      infile <- gsub(".gz", "", infile)
    }       
  }
  talk("Run pauda")
  ## get the paths to programs
  pauda_in_file <- ifelse(par_list["paired"],
                          gsub(".blastx", "_combinded_mod.fa", outfile), infile)
  pauda_out_file <- gsub(".blastx", "_pep.blastx", outfile)
  pauda_log_file <- file.path(par_list["tmp_dir"], "pauda.log")
  ## start default pauda
  pauda_cmd <- paste(program_list["pauda"],
                     pauda_in_file,
                     pauda_out_file,
                     par_list["amino_index"],
                     "1>", pauda_log_file, "2>&1"
                     )
  runCMD(pauda_cmd)
  ## get count df
  talk("Process pauda hits")
  pauda_hits_df <- get_pauda_hits(file = gsub(".blastx", "_pep.blastx", outfile),
                                  par_list) 
  ## clean everything
  unlink(dir0(dirname(outfile), "combinded"))
  ## count and return
  if("tbl_df" %in% class(pauda_hits_df)) {
    hits_table <- table(pauda_hits_df$genebank_id)
    hits_df <- arrange(tibble(genebank_id = names(hits_table), count_prot = hits_table),
                       desc(count_prot))
    hits_df <- filter(hits_df, count_prot > 5)
    return(list(count_df = hits_df,
                alignment_df = pauda_hits_df))
  } else {
    return(empty_map_pep_list())
  }
}
