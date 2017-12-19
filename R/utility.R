##' Test
##'
##' Test
##' @title Test 
##' @param batch_fastq_files 
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

##' Test
##'
##' Test
##' @title Test 
##' @param mat 
##' @param sep 
##' @return matrix
##' @author Jochen Kruppa
##' @export
collapse_df <- function(mat, sep = "_"){
  apply(mat, 1, function(x) str_c(x, collapse = sep))
}

##' Test
##'
##' Test
##' @title TEst
##' @param sample 
##' @param batch 
##' @param in_dir 
##' @param out_dir 
##' @param name 
##' @param gz 
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

##' Test
##'
##' Test
##' @title Test 
##' @param file 
##' @return logical
##' @author Jochen Kruppa
##' @export
is.named <- function(file){
    if(is.null(names(file))) {stop("Named file is needed")}
}

##' Test
##'
##' TEst
##' @title Test 
##' @param fqFile 
##' @param faFile 
##' @return NULL
##' @author Jochen Kruppa
##' @export
fq2fa <- function(fqFile, faFile) {
  runCMD(paste(program_list["seqtk"], "seq -A", fqFile, ">", faFile))
}


##' Test
##'
##' Test
##' @title Test 
##' @param dir 
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

##' Test
##'
##' Test
##' @title Test
##' @param map_dna_list 
##' @param map_pep_list 
##' @return data.frame
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

##' Test
##'
##' Test
##' @title Test 
##' @param sample_in 
##' @param out_file 
##' @param par_list 
##' @param proc_start_tm 
##' @return data.frame
##' @author Jochen Kruppa
##' @export
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
  if(length(names(sample_in)) >= 1){
    sample_name <- basename(out_file)
  }
  sample_info <- tibble(sample_name = sample_name,
                        num_reads = num_reads,
                        sequencing = ifelse(par_list["paired"], "paired-end", "single-end"),
                        host,
                        num_qc_reads = num_qc_reads)
  sample_info$run_time <- format(.POSIXct((proc.time() - proc_start_tm)[3], tz = "GMT"), "%H:%M:%S")
  return(sample_info)
}

##' Test
##'
##' Test
##' @title Test 
##' @param file 
##' @param force 
##' @return data.frame
##' @author Jochen Kruppa
##' @export
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
