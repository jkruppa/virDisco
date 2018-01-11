##' Small wrapper to extract the time from Sys.time()
##'
##' Small wrapper to extract the time from Sys.time()
##' @title Small wrapper to extract the time from Sys.time() 
##' @param ... No parameters are needed
##' @return System time in %b %d %X
##' @author Jochen Kruppa
TIME <- function(...) format(Sys.time(), "%b %d %X")


##' Advanced message() functionality by time stamp
##'
##' Advanced message() functionality by time stamp using the function
##' \code{\link{TIME}}.
##' @title Advanced message functionality by time stamp
##' @param ... Character string
##' @return Message with time stamp
##' @author Jochen Kruppa
##' @export
talk <- function(...) {
  require(stringr)
  message(str_c(TIME(), " ..... "), ...)
}

##' The function divides a DNAStringSet into a list of smaller objects given a threshold.
##'
##' Sometimes programs can not handle to long sequences. The functions
##' allows to build up chunks without cutting sequences. If the
##' threshold is exceeded by adding a sequence, the exceeding sequence
##' will be put in the next chunk.
##' @title Divide sequences into chunks
##' @param seqs DNAStringSet object with sequences
##' @param threshold maximal summed width of sequences in a chunk
##' @return list() of chunks with the width less the threshold
##' @author Jochen Kruppa
##' @export
get_seq_chunks_by_threshold <- function(seqs, threshold){
    ## get the file size
    seqList <- list()
    startCnt <- 1
    endCnt <- 1
    writeCnt <- 1
    repeat{
        tmpSize <- sum(width(seqs[startCnt:endCnt]))
        if(tmpSize >= threshold){
            endCnt <- endCnt - 1
            seqList[[writeCnt]] <- seqs[startCnt:endCnt]
            writeCnt <- writeCnt + 1
            startCnt <- endCnt + 1
        }
        if(endCnt >= length(seqs)) {
           seqList[[writeCnt]] <- seqs[startCnt:endCnt]
           return(seqList)
           break
        }
        endCnt <- endCnt + 1
    }
    names(seqList) <- str_pad(seq_along(seqList), 3, pad = "0")
    return(seqList)
}


##' This an experimental function, which is very specific
##'
##' In the daily workflow, some descriptions must be generated. This
##' is one of the functions.
##' @title Plot the processed files in a dir
##' @param batch_fastq_files A list of files to be parsed
##' @return ggplot object
##' @author Jochen Kruppa
##' @export
plot_time_sample <- function(batch_fastq_files) {
  require(scales)
  batch_raw_dates <- gsub("batch_", "", str_split(batch_fastq_files, "/", simplify = TRUE)[,6])
  batch_dates_df <- na.omit(tibble(dates = as.Date(batch_raw_dates, "%Y%m%d")))
  n_count <- nrow(batch_dates_df)/2
  p <- ggplot(batch_dates_df, aes(x = dates)) + geom_bar(aes(y = ..count../2)) +
    scale_x_date(labels = date_format("%b-%Y"),
                 breaks = unique(batch_dates_df)$dates) +
    annotate("text", x = unique(batch_dates_df)$dates[1], y = 25,
             label = str_c("Run samples: ", n_count), hjust = 0) +
    theme_bw() + ylab("Count")
  p <- p + theme(axis.text.x = element_text(angle = 60, hjust = 1))
  return(p)
}

##' Helper function for collapsing a matrix
##'
##' Helper function for collapsing a matrix
##' @title Collapse a matrix by row 
##' @param mat A matrix
##' @param sep A seperator for collapsing
##' @return matrix
##' @author Jochen Kruppa
##' @export
collapse_df <- function(mat, sep = "_"){
  apply(mat, 1, function(x) str_c(x, collapse = sep))
}

##' The function gets the file names and locations of a single sample
##'
##' The package virDisco needs a very specific file
##' declaration. Therefore these functions produces the structure.
##' @title Gather the files for one sample run
##' @param sample Name of the sample, which is also part of the file
##'   name
##' @param batch Name of the batch [default = NULL]
##' @param in_dir Name of the directory, where the fastq files are
##'   stored
##' @param out_dir Name of the directory, where the files should be
##'   written
##' @param name Name of the sample [expired]
##' @param gz Are the fastq files packed or not? [default = FALSE]
##' @return file.list
##' @author Jochen Kruppa
##' @export
get_sub_files <- function(sample, batch = NULL, in_dir = NULL, out_dir = NULL, name = NULL,
                          gz = FALSE) {
  if(is.null(name)){
    sample_name <- sample
  } else {
    sample_name <- name
  }
  ## in files
  if(is.null(batch) & !is.null(in_dir)){
    sample_in_dir <- dir0(in_dir)
  } else {
    sample_in_dir <- dir0(file.path(dataDir, batch, "input"))
  }
  if(any(file_ext(sample_in_dir) != "gz" & gz)){
    talk("Found unzipped input files, run gzip")
    l_ply(sample_in_dir, function(x) runCMD(paste("gzip", x)))
  }
  sample_in <- setNames(list(grep(str_c(sample, ".*.fastq.gz"), sample_in_dir,
                                         value = TRUE)), sample_name)
  if(length(unlist(sample_in)) == 1) {
    sample_in <- llply(sample_in, function(x) {names(x) <-  c("R1"); x})    
  } else {
    sample_in <- llply(sample_in, function(x) {names(x) <-  c("R1", "R2"); x})
  }
  ## out files
  if(is.null(out_dir)) {
    sample_proc_dir <- file.path(dataDir, batch, "proc")
  } else {
    sample_proc_dir <- out_dir
  }
  sample_paired_proc <- setNames(file.path(sample_proc_dir, sample_name, sample_name),
                                 sample_name)
  if(!file.exists(dirname(sample_paired_proc))) dir.create(dirname(sample_paired_proc))
  return(list(in_file = sample_in, out = sample_paired_proc))
}

##' Small tester if a object or file is named
##'
##' Small tester if a object or file is named
##' @title Small tester if a object is named
##' @param file A file string
##' @return logical
##' @author Jochen Kruppa
##' @export
is.named <- function(file){
    if(is.null(names(file))) {stop("Named file is needed")}
}

##' A wrapper for seqtk to convert fastq to fasta
##'
##' The external program seqtk must be installed and the path must be
##' set in the program_list().
##' @title A wrapper for seqtk to convert fastq to fasta
##' @param fqFile File path to fastq infile
##' @param faFile File path to fasta outfile
##' @return NULL
##' @author Jochen Kruppa
##' @export
fq2fa <- function(fqFile, faFile) {
  runCMD(paste(program_list["seqtk"], "seq -A", fqFile, ">", faFile))
}


##' A [expired] functionality to get the files from a directonary
##'
##' The function is _expired_. Use the function
##' \code{\link{get_sub_files}} instead.
##' @title Get all the files prepared in a directory
##' @param dir Directonary to be read
##' @return file.list
##' @author Jochen Kruppa
##' @export
prepare_file_list <- function(dir){
  raw_fq_paired <- grep(".fastq.gz", dir0(dir), value = TRUE)
  id_raw <- unique(str_split(basename(raw_fq_paired), "_", simplify = TRUE)[,1])
  if(any(grepl("Undetermined", id_raw))) {
    id <- id_raw[-grep("Undetermined", id_raw)]
  } else {
    id <- id_raw
  }
  ## name the files
  fq_paired <- llply(id, function(x){
    files <- grep(x, raw_fq_paired, value = TRUE)
    names(files) <- gsub(".*_(R\\d+)_.*", "\\1", files)
    return(files)
  })
  names(fq_paired) <- id
  return(fq_paired)
}

##' The hits from the DNA and AA mapping are merged and ordered by their combined rank.
##'
##' The hits from the DNA and AA mapping are merged and ordered by their combined rank.
##' @title Order the hits of the mapping
##' @param map_dna_list List with the DNA read hits from \code{\link{map_dna_ref}}
##' @param map_pep_list List with the AA read hits from \code{\link{map_pep_ref}}
##' @return Ordered data.frame with the counts
##' @author Jochen Kruppa
##' @export
ord_findings <- function(map_dna_list, map_pep_list) {
  dna_count_df <- arrange(map_dna_list$count_df, desc(count_read))
  dna_count_df$rank_dna <- 1:nrow(dna_count_df)
  pep_count_df <- arrange(map_pep_list$count_df, desc(count_prot))
  pep_count_df$rank_pep <- 1:nrow(pep_count_df)
  pep_count_df$genebank_id <- gsub(">", "", pep_count_df$genebank_id) 
  count_df <- left_join(dna_count_df, pep_count_df, by = c("genebank_id" = "genebank_id"))
  count_ord_df <- arrange(count_df, rank_dna + rank_pep) 
  return(count_ord_df)
}

##' A function to collect the sample information from the pipeline run.
##'
##' This is a internal function, which collects all information on the sample after the pipeline has run. 
##' @title Internal function to extract the sample information
##' @param sample_in Infile of the sample, see \code{\link{get_sub_files}}.
##' @param out_file Outfile of the sample, see \code{\link{get_sub_files}}.
##' @param par_list Parameter list, see \code{\link{set_par_list}}.
##' @param proc_start_tm String of the start time of the pipeline [internal]
##' @return data.frame
##' @author Jochen Kruppa
build_sample_info <- function(sample_in, out_file = NULL, par_list, proc_start_tm){
  talk("[SAMPLE INFO] Build sample info")
  talk("[SAMPLE INFO] Get the number of raw reads")
  if(par_list["qc"]) {
    trim_log_file <- dir0(par_list["tmp_dir"], str_c(basename(out_file), "_trim.log"))
    trim_info_df <- read.table(trim_log_file, header = TRUE)
    num_reads <- trim_info_df$num_reads
    num_qc_reads <- trim_info_df$num_qc_reads
  } else {
    ## if gz change to zcat
    cat <- ifelse(any(grepl(".gz$", unlist(sample_in))), "zcat", "cat")
    ## get the reads
    num_reads <- try(system(paste("echo $(", cat,
                                  unlist(sample_in)[1],
                                  "| wc -l)/4 | bc"),
                            intern = TRUE))
    ##
    num_qc_reads <- "No quality check and trimming executed"
  }
  if(par_list["check_host"]){
    talk("[SAMPLE INFO] Check host")
    host <- check_host(sample_in)
  } else {
    host <- "no host checked"    
  }
  sample_info <- tibble(sample_name = names(sample_in),
                        num_reads = num_reads,
                        sequencing = ifelse(par_list["paired"], "paired-end", "single-end"),
                        host,
                        num_qc_reads = num_qc_reads)
  sample_info$run_time <- format(.POSIXct((proc.time() - proc_start_tm)[3], tz = "GMT"), "%H:%M:%S")
  return(sample_info)
}

##' A experimental function, which checks the host using blastn
##'
##' This very experimental function uses blast and a subset of 100
##' read to determine the host of the sample. This works sometimes,
##' but it is very time consuming, with sometimes little gain.
##' @title Experimental function to check the host by blastn
##' @param file Infile of the sample, see \code{\link{get_sub_files}}.
##' @param force Should all files be overwritten?
##' @return data.frame
##' @author Jochen Kruppa
check_host <- function(file, force = FALSE){
  talk("Check host by blasting random reads (n = 100)")
  ## fles for blast
  blast_sub_fq <- file.path(tmpDir, gsub(".fastq|.fq.*", "_sub.fq", basename(unlist(file)[1])))
  blast_sub_fa <- gsub("fq.*", "fa", blast_sub_fq)
  blast_sub_out <- gsub("fq.*", "blast.out", blast_sub_fq)
  ## only run blast if the file does not exists
  if(!file.exists(blast_sub_out) || force){
    ## draw 100 reads from a fastq file
    runCMD(paste(program_list["seqtk"], "sample", unlist(file)[1], "100 >", blast_sub_fq))
    ## convert fq to fa
    fq2fa(blast_sub_fq, blast_sub_fa)
    ## run parallel blast
    multcoreBlastCMD <- paste("cat", blast_sub_fa, "| parallel --block 1k --recstart '>'",
                              "--pipe", program_list["blastn"],
                              "-db nt",
                              "-num_descriptions 1",
                              "-num_alignments 0", 
                              "-query -", ## dash is important!
                              ">", blast_sub_out)
    runCMD(multcoreBlastCMD)
  }
  ## get the blast results                   
  talk("Process blast output lines")
  rawBlastOut <- readLines(blast_sub_out)
  hitLinesPos <- grep("Sequences producing significant alignments", rawBlastOut) + 2
  hitLines <- unlist(llply(strsplit(rawBlastOut[hitLinesPos], "  "), function(x) x[2]))
  ## split the lines and get the most common
  split <- str_split(hitLines, " ", simplify = TRUE)[,1:2]
  hit_table <- table(str_c(split[,1], split[,2], sep = " "))
  hit_max <- sort(hit_table, decreasing = TRUE)[1:2] 
  return(str_c(names(hit_max[1]), " (", hit_max[1], "%)",
               " or ",
               names(hit_max[2]), " (", hit_max[2], "%)"))
}
