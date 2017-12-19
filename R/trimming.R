##' Test
##'
##' Test
##' @title Test 
##' @param inFile 
##' @param outFile 
##' @param leading 
##' @param trailing 
##' @param minlength 
##' @param illuminaclip 
##' @param log_file 
##' @return Test
##' @author Jochen Kruppa
##' @export
fastq_trimmer <- function (inFile, outFile, leading = 10, trailing = 10, minlength = 50, 
                           illuminaclip, log_file) 
    {
      if(length(inFile) == 2) {
        paired_outFile <- gsub("trimmed", "trimmed.paired", outFile)
        unpaired_outFile <- gsub("trimmed", "trimmed.unpaired", outFile)
        fastq_trimmer_CMD <- paste("java -jar", program_list["trimmomatic"], 
                                   "PE",
                                   "-threads", par$nCores,
                                   "-phred33",
                                   inFile["R1"], 
                                   inFile["R2"],
                                   paired_outFile["R1"],
                                   unpaired_outFile["R1"], 
                                   paired_outFile["R2"],
                                   unpaired_outFile["R2"],
                                   str_c("ILLUMINACLIP:", illuminaclip, ":2:30:10"),
                                   str_c("LEADING:", leading),
                                   str_c("TRAILING:", trailing),
                                   str_c("MINLEN:", minlength),
                                   "SLIDINGWINDOW:4:15",
                                   "2>", log_file)
    }
    else {
      fastq_trimmer_CMD <- paste("java -jar", program_list["trimmomatic"], 
                                 "SE",
                                 "-threads", par$nCores,
                                 "-phred33",
                                 inFile, 
                                 outFile,
                                 str_c("ILLUMINACLIP:", illuminaclip, ":2:30:10"), 
                                 str_c("LEADING:", leading),
                                 str_c("TRAILING:", trailing), 
                                 str_c("MINLEN:", minlength),
                                 "SLIDINGWINDOW:4:15",
                                 "2>", log_file)
    }
    runCMD(fastq_trimmer_CMD)
}


##' Main fastq quality function
##'
##' Main fastq quality function. Calls the function
##' fastq_quality_filter and fastq_trimmer. Stores files in a tempDir.
##' @title Main fastq quality function
##' @param inFile Fastq infile
##' @param illumninaclip Remove the adapter from illumina by TruSeq2
##'   or TruSeq3 protocol
##' @param par_list 
##' @param leading Numeric. Remove leading low quality or N bases
##'   (below quality 10 [default])
##' @param trailing Numeric. Remove trailing low quality or N bases
##'   (below quality 10 [default])
##' @param minlength Min read length to keep. 50 [default]
##' @param log_file write the trimmomatic output to file
##' @return NULL
##' @author Jochen Kruppa
##' @export
fastq_quality_control <- function(inFile,
                                  illumninaclip,
                                  par_list,
                                  leading = 10, 
                                  trailing = 10,
                                  minlength = 50)
{
  ## talk("Write all temporal files to /home/temp")
  log_file <- file.path(par_list["tmp_dir"],
                        str_c(names(inFile), "_trim.log"))
  tmp_in_file <- unlist(inFile)
  names(tmp_in_file) <- gsub(".*\\.", "", names(tmp_in_file))
  if(length(tmp_in_file) == 2){
    tmpTrimFq <- file.path(par_list["tmp_dir"], gsub("fastq|fq", "trimmed.fq", basename(tmp_in_file)))
    names(tmpTrimFq) <- names(tmp_in_file)
    talk("[TRIMMING] Start trimming")
    fastq_trimmer(inFile = tmp_in_file, outFile = tmpTrimFq, leading, trailing, minlength, 
                  illumninaclip, log_file)
    out_file <- list(setNames(laply(tmpTrimFq, function(x) gsub("trimmed.", "trimmed.paired.", x)),
                              c("R1", "R2")))
    names(out_file) <- names(inFile)
    ## min length is really 50
    talk(str_c("[TRIMMING] Remove reads shorter than ", par_list["min_num_reads"], "bp"))
    trim_fq_list <- llply(unlist(out_file), readFastq)
    num_reads <- sum(laply(trim_fq_list, length))
    long_pos <- intersect(which(width(trim_fq_list[[1]]) >= par_list["min_num_reads"]),
                          which(width(trim_fq_list[[2]]) >= par_list["min_num_reads"]))
    num_qc_reads <- 2 * length(long_pos)
    mean_read_length <- mean(width(trim_fq_list[[1]][long_pos]))
    unlink(unlist(out_file))
    l_ply(seq_along(trim_fq_list), function(i) writeFastq(trim_fq_list[[i]][long_pos],
                                                          unlist(out_file)[i]),
          .parallel = TRUE)
    ## write trim_log
    trim_log <- tibble(num_reads, num_qc_reads, mean_read_length)
    write_delim(trim_log, log_file)    
  } else {
    ## this is redundant I know, but single reads are not fully tested yet
    tmpTrimFq <- file.path(par_list["tmp_dir"], gsub("fastq|fq", "trimmed.fq", basename(tmp_in_file)))
    names(tmpTrimFq) <- names(tmp_in_file)
    talk("[TRIMMING] Start trimming")
    fastq_trimmer(inFile = tmp_in_file, outFile = tmpTrimFq, leading, trailing, minlength, 
                  illumninaclip, log_file)     
    out_file <- list(setNames(laply(tmpTrimFq, function(x) gsub("trimmed.", "trimmed.", x)),
                              c("R1")))
    names(out_file) <- names(inFile)
    ## min length is really 50
    talk(str_c("[TRIMMING] Remove reads shorter than ", par_list["min_num_reads"], "bp"))
    trim_fq_list <- llply(unlist(out_file), readFastq)
    num_reads <- sum(laply(trim_fq_list, length))
    long_pos <- which(width(trim_fq_list[[1]]) >= par_list["min_num_reads"])
    num_qc_reads <- length(long_pos)
    mean_read_length <- mean(width(trim_fq_list[[1]][long_pos]))
    unlink(unlist(out_file))
    l_ply(seq_along(trim_fq_list), function(i) writeFastq(trim_fq_list[[i]][long_pos],
                                                          unlist(out_file)[i]),
          .parallel = TRUE)
    ## write trim_log
    trim_log <- tibble(num_reads, num_qc_reads, mean_read_length)
    write_delim(trim_log, log_file)    
  }
  return(out_file)
}
