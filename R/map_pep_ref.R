##' This function calls PANDAseq and PAUDA for mapping the translated DNA reads
##'
##' The function is the mapping control function for the AA
##' mapping. Most importantly, paired reads are merged by PAUDAseq
##' first. The mapping and translation of the sequence reads to amino
##' reads is done by PAUDA. Finally, the blastx output is parsed by
##' \code{\link{get_pauda_hits}}.
##' @title Amino mapping function
##' @param infile File path to the fastq file, better controlled by
##'   parameter given by the par_list(), see
##'   \code{\link{set_par_list}}
##' @param outfile File path to the results dir
##' @param par_list Parameter given by the par_list(), see
##'   \code{\link{set_par_list}}
##' @return pep_alignment_df list
##' @author Jochen Kruppa
##' @export
##' @examples
##' data(NHS_10001_map_aa_list)
map_pep_ref <- function(infile, outfile, par_list){
  ## mapping to pep reference
  if(par_list["paired"]) {
    infile <- unlist(infile)
    talk("[AMINO MAPPING] Run pandaseq to build up one fastq file")
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
  talk("[AMINO MAPPING] Run pauda")
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
  talk("[AMINO MAPPING] Process pauda hits")
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
