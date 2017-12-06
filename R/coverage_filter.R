##' Desc
##'
##' Detail
##' @title Test 
##' @param hit_id 
##' @param map_dna_list 
##' @param par_list 
##' @param out_file 
##' @return tibble
##' @author Jochen Kruppa
##' @export
coverage_filter <- function(hit_id, map_dna_list, par_list, out_file){
  viral_length <- collect(filter(par_list["species_info"], genebank_id %in% hit_id)) %>%
    select(genebank_id, length = Length)
  coverage_tbl <- ldply(hit_id, function(x) {
    viral_width <- filter(map_dna_list$alignment_df, genebank_id == x) %$%
      IRanges::reduce(IRanges(start = .$pos_start, end = .$pos_end)) %>%
      width(.)
    viral_read_length_mapq <- select(filter(map_dna_list$alignment_df, genebank_id == x),
                                     genebank_id, read_length, mapq)
    coverage <- sum(viral_width)/filter(viral_length, genebank_id == x)$length
    tibble(genebank_id = x,
           coverage = round(coverage, 4),
           mean_mapq = mean(viral_read_length_mapq$mapq),
           mean_read_length = mean(viral_read_length_mapq$read_length),
           hits = length(viral_width))
  }, .parallel = TRUE)
  png(str_c(out_file, "_coverage.png"), width = 18, height = 18, res = 300, units = "cm")
  p <- ggplot(coverage_tbl, aes(coverage)) +
    geom_density(fill = par$cbbPalette[3], alpha = 0.6) +
    theme_bw() +
    xlab("Coverage of the reference genome") + ylab("Density") +
    ggtitle(str_c(basename(out_file), " - Coverage"))  +
    geom_vline(xintercept = 0.05, color = par$cbbPalette[7])
  print(p)
  dev.off()
  png(str_c(out_file, "_mapq_hist.png"), width = 18, height = 18, res = 300, units = "cm")
  p <-   ggplot(coverage_tbl, aes(mean_mapq)) +
    geom_density(fill = par$cbbPalette[3], alpha = 0.6) +
    theme_bw() + xlim(c(0, 255)) +
    xlab("Mapping quality") + ylab("Density") +
    ggtitle(str_c(basename(out_file), " - Mean mapping quality"))    
  print(p)
  dev.off()
  png(str_c(out_file, "_read_length_hist.png"), width = 18, height = 18, res = 300, units = "cm")
  dens <- density(coverage_tbl$mean_read_length)
  max_y <- which.max(dens$y)
  x_break <- c(dens$x[max_y], round(seq(from = min(dens$x), to = max(dens$x), length.out = 5)))
  p <- ggplot(coverage_tbl, aes(mean_read_length)) +
    geom_density(fill = par$cbbPalette[3], alpha = 0.6) +
    theme_bw() + xlab("Read length") + ylab("Density") +
    geom_vline(xintercept = dens$x[max_y], color = par$cbbPalette[7]) +
    ggtitle(str_c(basename(out_file), " - Mean read length")) +
    scale_x_continuous(breaks = x_break)
  print(p)
  dev.off()  
  return(coverage_tbl)
}

