##' Test
##'
##' Test
##' @title Test
##' @param infile 
##' @param outfile 
##' @param par_list 
##' @param min_hit 
##' @param all 
##' @return data.frame
##' @author Jochen Kruppa
##' @export
map_dna_ref <- function(infile, outfile, par_list, min_hit = 5, all = FALSE){
  ## mapping to dna reference
  bam_file <- gsub(".sam", "_dna.bam", outfile)
  if(!file.exists(bam_file)){
    switch(par_list["mapper"],
           "bowtie" = {
             bowtie2CMD <- get_bowtie2_cmd(inFile = infile,
                                           samOutFile = gsub(".sam", "_dna.sam", outfile),
                                           referenceDir = par_list["genome_index"],
                                           method = "end-to-end",
                                           p = par_list["ncore"],
                                           all)
             runCMD(bowtie2CMD)
           },
           "star" = {
             l_ply(unlist(infile), function(x) {
               runCMD(paste("gunzip -k -f", x))
             })
             star_CMD <- paste(STAR,
                               "--runThreadN", par_list["ncore"],
                               "--genomeDir", par_list["genome_index"],
                               "--readFilesIn", str_c(gsub(".gz", "", unlist(infile)), collapse = " "),
                               "--outFileNamePrefix", gsub(".sam", "_", outfile),
                               "1>", file.path(tmpDir, "star_mapper.log"), "2>&1")
             runCMD(star_CMD)
             ## rename
             file.rename(str_replace(outfile, ".sam", "_Aligned.out.sam"),
                         str_replace(outfile, ".sam", "_dna.sam"))
             ##
             unlink(gsub(".gz", "", unlist(infile)))
             unlink(list.files(dirname(infile[[1]][1]), "STARtmp", full.names = TRUE), recursive = TRUE)
           })
    ## convert sam to bam  
    sam2bam_serial(sam_files = gsub(".sam", "_dna.sam", outfile), keep = FALSE)
  }
  ## read in bam files
  alignment_bam <- Reduce(c, scanBam(bam_file))
  ## get read / reference information (qname is read, rname is reference)
  read_ref_tbl <- tibble(read = alignment_bam$qname, ref = alignment_bam$rname) %>%
    mutate(read = str_replace(read, "R1_|R2_", "")) %>%
    mutate(read = ifelse(str_detect(read, "decoy"), "decoy", "sample")) %>%
    mutate(ref = ifelse(str_detect(ref, "decoy"), "decoy", "viral")) %>%
    mutate(match = str_c(read, "_", ref))
  ## build up the 2x2 table
  read_ref_table <- table(reads = read_ref_tbl$read, ref = read_ref_tbl$ref)
  multi_map_rate <- sum(read_ref_table["decoy", ]) / (2 * par_list["num_decoy_reads"])
  true_positive <- read_ref_table["decoy", "decoy"] / sum(read_ref_table["decoy", ])
  false_positive <- (read_ref_table["decoy", "viral"] + read_ref_table["sample", "decoy"]) / sum(read_ref_table)
  ## get the aligment data frame for plotting
  qname_raw <- str_split(alignment_bam$qname, ":", simplify = TRUE)
  alignment_df <- tibble(genebank_id = alignment_bam$rname,
                             qname = qname_raw[, ncol(qname_raw)],
                             strand = alignment_bam$strand,
                             read_length = alignment_bam$qwidth,
                             pos_start = alignment_bam$pos,
                             pos_end = pos_start + read_length - 1,
                             mapq = alignment_bam$mapq,
                             seq_id = str_c(genebank_id, qname, pos_start, pos_end, sep = "_"))
  ## get the raw read seqs
  seqs_raw <- alignment_bam$seq
  names(seqs_raw) <- with(alignment_raw_df, str_c(genebank_id, qname, pos_start, pos_end, sep = "_"))
  ## filter the seqs
  seqs <- seqs_raw[names(seqs_raw) %in% alignment_df$seq_id]
  seq_file <- gsub(".sam", "_dna_seqs.RDS", outfile)
  saveRDS(seqs, seq_file)
  ## get the counts
  alignment_table <- table(alignment_df$genebank_id) 
  alignment_table_sorted <- sort(alignment_table, decreasing = TRUE)
  alignment_table_clean <- alignment_table_sorted[alignment_table_sorted > min_hit]
  ## get the count data frame
  count_df <- setNames(ldply(alignment_table_clean), c("genebank_id", "read_count")) %>%
    extract(-str_which(.$genebank_id, "decoy"),)
  ## return list with counts and aligment information
  return(list(count_df = select(count_df, genebank_id, count_read = read_count),
              alignment_df = alignment_df,
              stats = c(multi_map_rate = multi_map_rate,
                        true_positive = true_positive,
                        false_positive = false_positive)))
}
