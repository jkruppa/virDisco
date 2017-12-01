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
##' @return NULL
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
##' @param tmpDir Where to store the temporal files. Files will be
##'   removed automatically
##' @param leading Numeric. Remove leading low quality or N bases
##'   (below quality 10 [default])
##' @param trailing Numeric. Remove trailing low quality or N bases
##'   (below quality 10 [default])
##' @param minlength Min read length to keep. 50 [default]
##' @param illumninaclip Remove the adapter from illumina by TruSeq2
##'   or TruSeq3 protocol
##' @param log_file write the trimmomatic output to file
##' @param outFile Fastq outfile
##' @return NULL
##' @author Jochen Kruppa
##' @export
fastq_quality_control <- function(inFile, tmpDir, leading = 10, 
                                  trailing = 10, minlength = 50,
                                  illumninaclip, log_file)
{
  ## talk("Write all temporal files to /home/temp")
  tmp_in_file <- unlist(inFile)
  names(tmp_in_file) <- gsub(".*\\.", "", names(tmp_in_file))
  if(length(tmp_in_file) == 2){
    tmpTrimFq <- file.path(tmpDir, gsub("fastq|fq", "trimmed.fq", basename(tmp_in_file)))
    names(tmpTrimFq) <- names(tmp_in_file)
    talk("Start trimming")
    fastq_trimmer(inFile = tmp_in_file, outFile = tmpTrimFq, leading, trailing, minlength, 
                  illumninaclip, log_file)
    out_file <- list(setNames(laply(tmpTrimFq, function(x) gsub("trimmed.", "trimmed.paired.", x)),
                              c("R1", "R2")))
    names(out_file) <- names(inFile)
  } else {
    ## this is redundant I know, but single reads are not fully tested yet
    tmpTrimFq <- file.path(tmpDir, gsub("fastq|fq", "trimmed.fq", basename(tmp_in_file)))
    names(tmpTrimFq) <- names(tmp_in_file)
    talk("Start trimming")
    fastq_trimmer(inFile = tmp_in_file, outFile = tmpTrimFq, leading, trailing, minlength, 
                  illumninaclip, log_file)     
    out_file <- list(setNames(laply(tmpTrimFq, function(x) gsub("trimmed.", "trimmed.", x)),
                              c("R1")))
    names(out_file) <- names(inFile)
  }
  return(out_file)
}
