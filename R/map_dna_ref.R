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
                                           p = 32,
                                           all)
             runCMD(bowtie2CMD)
           },
           "star" = {
             l_ply(unlist(infile), function(x) {
               runCMD(paste("gunzip -k -f", x))
             })
             star_CMD <- paste(STAR,
                               "--runThreadN", 32,
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
  ## get the aligment data frame for plotting
  qname_raw <- str_split(alignment_bam$qname, ":", simplify = TRUE)
  alignment_raw_df <- tibble(genebank_id = alignment_bam$rname,
                             qname = qname_raw[, ncol(qname_raw)],
                             strand = alignment_bam$strand,
                             read_length = alignment_bam$qwidth,
                             pos_start = alignment_bam$pos,
                             pos_end = pos_start + read_length - 1,
                             mapq = alignment_bam$mapq,
                             seq_id = str_c(genebank_id, qname, pos_start, pos_end, sep = "_"))
  alignment_df <- filter(alignment_raw_df, mapq > 0 & read_length > (0.25 * max(read_length)))
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
  count_df <- setNames(ldply(alignment_table_clean), c("genebank_id", "read_count"))
    ## get in the decoy information
  decoy_pos <- str_detect(count_df$genebank_id, "decoy")
  false_positive <- count_df[decoy_pos,]$read_count / sum(count_df[!decoy_pos,]$read_count)
  ## adjusted count by false positives
  count_df %<>% 
    mutate(read_count_adj = round(read_count * false_positive)) %>%
    extract(!decoy_pos,)
  ## return list with counts and aligment information
  return(list(count_df = select(count_df, genebank_id, count_read = read_count),
              alignment_df = alignment_df,
              false_positive = false_positive))
}
