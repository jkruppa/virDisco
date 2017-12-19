##' Test
##'
##' Test
##' @title Test 
##' @param genebank_id 
##' @param consensus_df 
##' @param dna_alignment_df 
##' @param pep_alignment_df 
##' @param aa_info_df 
##' @param viral_info_df 
##' @param sample_info 
##' @param pdf_file 
##' @param par_list 
##' @return NULL
##' @author Jochen Kruppa
##' @export
mapping_dna_plot <- function(genebank_id,
                             consensus_df = NULL,
                             dna_alignment_df,
                             pep_alignment_df,
                             aa_info_df,
                             viral_info_df,
                             sample_info = NULL,
                             pdf_file = NULL,
                             par_list)
{
  pep_alignment_df$ident_percent <- as.numeric(gsub(".*\\((.*)%\\)$", "\\1", pep_alignment_df$ident))/100
  dna_pep_hits <- union(unique(dna_alignment_df$genebank_id),
                        unique(filter(pep_alignment_df, ident_percent > 0.9)$genebank_id))
  dna_pep_hits <- dna_pep_hits[order(match(dna_pep_hits, genebank_id))][1:length(genebank_id)]
  ## create png_dir
  out_file <- gsub(".pdf", "", pdf_file)
  png_dir <- file.path(str_c(out_file, "_png"))
  if(!file.exists(png_dir)) dir.create(png_dir, showWarnings = FALSE)
  if(par_list["plot"]){
    png_files <- file.path(png_dir,
                           str_c(basename(out_file), "_",
                                 str_pad(0:length(dna_pep_hits), 2, pad = "0"), ".png"))
    png_file_df <- tibble(genebank_id = dna_pep_hits, png_file = png_files[-1])
  } else {
    png_files <- file.path(png_dir, str_c(basename(out_file), "_0", ".png"))
  }
  ## start png production
  if(!is.null(sample_info)){
    png(png_files[1], width = 18, height = 18, res = 300, units = "cm")
    p <- ggplot(tibble(x = 0:20, y = 0:20), aes(x, y)) + theme_bw() + geom_blank() +
      annotate("text", 0, 20, label = str_c("Sample: ", sample_info$sample_name),
               hjust = 0, size = 10) +
      annotate("text", 0, 18, label = str_c('paste(bold("Date: "), \"', date(), "\")"),
               hjust = 0, parse = TRUE) +
      annotate("text", 0, 17, label = str_c('paste(bold("Number of raw reads: "), ',
                                            sample_info$num_reads, ")"),
               hjust = 0, parse = TRUE) +
      annotate("text", 0, 16, label = str_c('paste(bold("Number of reads after QC: "), \"',
                                            sample_info$num_qc_reads, "\")"),
               hjust = 0, parse = TRUE) +
      annotate("text", 0, 15,
               label = str_c('paste(bold("Sequence type: "), ', sample_info$sequencing, ")"),
               hjust = 0, parse = TRUE) +
      annotate("text", 0, 14,
               label = str_c('paste(bold("Host species: "), \"', sample_info$host, "\")"),
               hjust = 0, parse = TRUE) +
      annotate("text", 0, 13,
               label = str_c('paste(bold("Multi map rate: "), \"',
                             round(par_list["map_dna_stats"]$multi_map_rate, 2), "\")"),
               hjust = 0, parse = TRUE) +
      annotate("text", 0, 12,
               label = str_c('paste(bold("True positive rate: "), \"',
                             round(par_list["map_dna_stats"]$true_positive, 2), "\")"),
               hjust = 0, parse = TRUE) +
      annotate("text", 0, 11,
               label = str_c('paste(bold("False positive rate: "), \"',
                             round(par_list["map_dna_stats"]$false_positive, 2), "\")"),
               hjust = 0, parse = TRUE) +
      annotate("text", 0, 10,
               label = str_c('paste(bold("Run time: "), \"', sample_info$run_time, "\")"),
               hjust = 0, parse = TRUE) +  
      ## bottom
      annotate("text", 20, 0,
               label = "The number of reads is estimated by trimmomatic",
               hjust = 1) +
      annotate("text", 20, 1,
               label = "The host species is determined from 100 randomly selected reads",
               hjust = 1) +
      theme(panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank()) 
    print(p)
    dev.off()
  }
  if(par_list["plot"]){
    l_ply(dna_pep_hits, function(x) {
    png(filter(png_file_df, genebank_id  == x)$png_file, width = 18, height = 18, res = 300, units = "cm")
    flag_dna_reduce <- FALSE
    flag_pep_reduce <- FALSE
    virus_info <- collect(filter(viral_info_df, genebank_id == x))
    tmp_df <- filter(dna_alignment_df, genebank_id == x)
    tmp_df$id <- seq(0.5, 8, length.out = nrow(tmp_df))
    ## get the gene positions
    tmp_gene_df <- collect(filter(aa_info_df, genebank_id == x))
    ## plot pep reads
    prot_start_map <- hashmap(tmp_gene_df$prot_id, as.numeric(tmp_gene_df$pos_start))
    tmp_pep_pos_df <- filter(pep_alignment_df, genebank_id == x)
    tmp_pep_df <- tibble(qname = tmp_pep_pos_df$qname,
                         pos_start = prot_start_map[[tmp_pep_pos_df$prot_id]] + 3 * tmp_pep_pos_df$pos_start,
                         pos_end = prot_start_map[[tmp_pep_pos_df$prot_id]] + 3 * tmp_pep_pos_df$pos_end)
    tmp_pep_df <- arrange(tmp_pep_df, pos_start) 
    tmp_pep_df$id <- seq(-0.5, -8, length.out = nrow(tmp_pep_df))
    ## num of reads
    num_dna <- nrow(tmp_df)
    num_pep <- nrow(tmp_pep_df)
    ## get the intersect of reads
    intersect_reads_num <- length(intersect(tmp_pep_df$qname, tmp_df$qname))
    ## get the overlap
    if(nrow(tmp_df) == 0 || nrow(tmp_pep_df) == 0){
      overlap_reads_num <- 0
    } else {
      pep_range <- IRanges(tmp_pep_df$pos_start, tmp_pep_df$pos_end)
      dna_range <- IRanges(tmp_df$pos_start, tmp_df$pos_end)
      if(nrow(tmp_df) >= 10000 || nrow(tmp_pep_df) >= 10000){
        pep_range <- IRanges::reduce(pep_range)
        dna_range <- IRanges::reduce(dna_range)
      }
      overlap_dna_reads_01 <- laply(seq_along(dna_range), function(i) {
        ifelse(countOverlaps(dna_range[i], pep_range) != 0, 1, 0)
      })
      overlap_pep_reads_01 <- laply(seq_along(pep_range), function(i) {
        ifelse(countOverlaps(pep_range[i], dna_range) != 0, 1, 0)
      })
      overlap_reads_num <- sum(overlap_dna_reads_01) + sum(overlap_pep_reads_01)
    }
    ## start plot
    p <- ggplot() 
    ## draw the genes
    if(nrow(tmp_gene_df) > 0){
      p <- p + geom_rect(data = tmp_gene_df,
                         aes(xmin = pos_start, xmax = pos_end,
                             ymin = -Inf, ymax = Inf),
                         fill = alpha("gray", 0.25)) +
        geom_segment(data = tmp_gene_df,
                     aes(x = pos_start, xend = pos_end,
                         y = 0, yend = 0),
                     color = "red", size = 2) 
    }
    ## draw the reads
    if(nrow(tmp_df) >= 5000){ 
      reduced_dna_range <- IRanges::reduce(IRanges(start = tmp_df$pos_start, end = tmp_df$pos_end))
      tmp_df <- tibble(pos_start = start(reduced_dna_range),
                       pos_end = end(reduced_dna_range),
                       id = seq(0.5, 8, length.out = length(reduced_dna_range)))
      num_dna_reduced <- length(reduced_dna_range)
      flag_dna_reduce <- TRUE
    }
    if(nrow(tmp_pep_df) >= 5000){
      reduced_pep_range <- IRanges::reduce(IRanges(start = tmp_pep_df$pos_start, end = tmp_pep_df$pos_end))
      tmp_pep_df <- tibble(pos_start = start(reduced_pep_range),
                           pos_end = end(reduced_pep_range),
                           id = seq(-0.5, -8, length.out = length(reduced_pep_range)))
      num_pep_reduced <- length(reduced_pep_range)
      flag_pep_reduce <- TRUE
    }  
    p <- p + geom_segment(data = tmp_df,
                          aes(x = pos_start, xend = pos_end, y = id, yend = id),
                          color = par$cbbPalette[6]) +
      ## draw the pep hits
      geom_segment(data = tmp_pep_df,
                   aes(x = pos_start, xend = pos_end, y = id, yend = id),
                   color = par$cbbPalette[4]) +
      ## draw the reference
      geom_segment(aes(x = 1,
                       xend = virus_info$Length,
                       y = 0, yend = 0)) +
      theme_bw()
    if(flag_dna_reduce){
      p <- p + annotate("text", x = 0, y = 10,
                        label = 'bold("To much reads to plot: Build and plot pseudo contigs")',
                        size = 4,
                        parse = TRUE, hjust = 0)      
      p <- p + annotate("text", x = 0, y = 9,
                        label = str_c('paste(bold("DNA"), " -- ', num_dna, " reads reduced to ",
                                      num_dna_reduced, " contigs\")"),
                        size = 4,
                        parse = TRUE, hjust = 0)
    } else {
      p <- p + annotate("text", x = 0, y = 10,
                        label = str_c('paste(bold("DNA"), " -- ', num_dna, " reads\")"),
                        size = 4,
                        parse = TRUE, hjust = 0)
    }
    if(nrow(tmp_gene_df) == 0){
      p <- p + annotate("text", x = 0, y = -10,
                        label = 'paste(bold("AMINO"), " -- No gene information provided by NCBI")',
                        size = 4,
                        parse = TRUE, hjust = 0)
    } else {
      if(flag_pep_reduce){
        p <- p + annotate("text", x = 0, y = -9,
                          label = 'bold("To much reads to plot: Build and plot pseudo contigs")',
                          size = 4,
                          parse = TRUE, hjust = 0)      
        p <- p + annotate("text", x = 0, y = -10,
                          label = str_c('paste(bold("AMINO"), " -- ', num_pep, " reads reduced to ",
                                        num_pep_reduced, " contigs\")"),
                          size = 4,
                          parse = TRUE, hjust = 0)
      } else {
        p <- p + annotate("text", x = 0, y = -10,
                          label = str_c('paste(bold("AMINO"), " -- ', num_pep, " reads\")"),
                          size = 4,
                          parse = TRUE, hjust = 0)
      }
    }
    p <- p +  xlab("Reference genome") + ylab("Single reads") +
      ggtitle(str_c(virus_info$genebank_id, " ", virus_info$Description),
              subtitle = str_c(num_dna + num_pep,
                               " DNA and AMINO reads map to the reference ",
                               "(", intersect_reads_num, " identical by name",
                               ", ", overlap_reads_num, " overlap by position)")) +
      scale_y_discrete(breaks = "")  +
      theme(legend.position = "none") #+
    ## get the consensus plot
    if(par_list["consensus"]) {
      tmp_con_df <- filter(consensus_df, genebank_id == x)
    } else {
      tmp_con_df <- tibble()
    }
    if(nrow(tmp_con_df) == 0) {
      par_list["consensus"] <- FALSE
    }
    if(par_list["consensus"]) {
      p_2 <- ggplot(data = tmp_con_df,
                    aes(x = pos, y = con_prob, colour = bool)) +
        geom_vline(aes(xintercept = pos, color = bool)) +
        ## geom_smooth(aes(group = 1), span = 0.7, color = "black", se = FALSE) +
        geom_path(aes(group = 1), color = "black", alpha = 0.5) +
        scale_color_manual(values = par$cbbPalette[c(7, 6)]) +
        scale_y_continuous(breaks = NULL) +
        xlab("") + ylab("Probability") +
        theme_bw()
      p_2 <- p_2 + theme(legend.position = "none")
    }
    ## both together
    if(par_list["consensus"]) {
      require(gridExtra)
      require(grid)
      suppressMessages(grid.arrange(p, p_2,
                                    ncol = 1, nrow = 2, widths = c(5), heights = c(4, 1)))
    } else {
      print(p)
    }
    dev.off()
    }, .progress = "text")
    ## put the png's into one pdf
  }
  talk("Transfrom png files to one pdf")
  png_files <- c(png_files[1],
                 str_c(out_file, "_coverage.png"),
                 str_c(out_file, "_read_length_hist.png"),
                 png_files[2:length(png_files)])
  png_plots <- llply(png_files, function(x) {
    rasterGrob(readPNG(x, native = FALSE), interpolate = FALSE)
  })
  pdf(pdf_file)
  l_ply(png_plots, function(x) print(ggplot() + annotation_custom(x)))
  dev.off()
}

