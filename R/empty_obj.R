##' Test
##'
##' Test
##' @title Test
##' @param ... 
##' @return data.frame
##' @author Jochen Kruppa
##' @export
empty_map_pep_list <- function(...){
  list(count_df = tibble(genebank_id = "",
                         count_prot = ""),
       alignment_df = tibble(qname = "",
                             prot_id = "",
                             pos_start = NA,
                             pos_end = NA,
                             genebank_id = NA,
                             ident = 0))
}
