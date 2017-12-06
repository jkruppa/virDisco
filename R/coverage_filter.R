##' Desc
##'
##' Detail
##' @title Test 
##' @param hit_id 
##' @param par_list 
##' @return tibble
##' @author Jochen Kruppa
##' @export
coverage_filter <- function(hit_id, par_list){
  viral_length <- collect(filter(par_list["species_info"], genebank_id %in% hit_id)) %>%
    select(genebank_id, length = Length)
  coverage_tbl <- ldply(hit_id, function(x) {
    viral_width <- filter(map_dna_list$alignment_df, genebank_id == x) %$%
      IRanges::reduce(IRanges(start = .$pos_start, end = .$pos_end)) %>%
      width(.) 
    coverage <- sum(viral_width)/filter(viral_length, genebank_id == x)$length
    tibble(genebank_id = x,
           coverage = round(coverage, 4),
           hits = length(viral_width))
  }, .parallel = TRUE)
  return(coverage_tbl)
}
