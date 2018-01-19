##' The function extracts the mapping information of PAUDA
##'
##' Internal function of \code{\link{map_pep_ref}} to extract the
##' information of the blastx output of the amino acid mapping by
##' PAUDA.
##' @title Internal parse function for PAUDA
##' @param file blastx input file from PAUDA output
##' @param par_list Parameter given by the par_list(), see
##'   \code{\link{set_par_list}}
##' @return data.frame
##' @author Jochen Kruppa
##' @export
get_pauda_hits <- function(file, par_list) {
  tmp_lines <- readLines(file)
  ## cut the queries into blocks
  query_start <- grep("^Query=", tmp_lines)
  query_end <- c(query_start[-1] - 2, length(tmp_lines))
  query_block_list <- llply(seq_along(query_start), function(i){
    tmp_lines[query_start[i]:query_end[i]]
  })  
  ## get the entries in a single block
  query_block_raw_df <- tbl_df(ldply(query_block_list, function(x){
    score_line <- grep(" Score", x, value = TRUE)
    ident_line <- grep(" Identities", x, value = TRUE)
    sbjct_line <- grep("Sbjct:", x, value = TRUE)
    query_line <- grep("Query:", x, value = TRUE)
    qname <- str_split(gsub("^Query= ", "", grep("^Query=", x, value = TRUE)), ":", simplify = TRUE)
    tibble(qname = qname[, length(qname)],
           prot_id = gsub("^>", "", str_split(grep("^>", x, value = TRUE),
                                              par_list["gen_prot_sep"], simplify = TRUE)[,1]),
           pos_start = as.numeric(gsub("Sbjct: (\\d+)\\s+.*", "\\1", sbjct_line)),
           pos_end = as.numeric(gsub(".*\\s+(\\d+)$", "\\1", sbjct_line)),
           score_bit = as.numeric(gsub(" Score = (.*) bits.*", "\\1", score_line)),
           score_expect = as.numeric(gsub(".*Expect = (.*)$", "\\1", score_line)),
           ident = gsub(" Identities = (.*), Pos.*", "\\1", ident_line),
           ident_pos = gsub(".* Positives = (.*), Gaps.*", "\\1", ident_line)
           ## aa_seq =  llply(gsub("Query: \\d+\\s+(\\w*)\\s\\d+", "\\1", query_line), AAString)
           )
  # }, .progress = "text"))
  }, .parallel = TRUE))
  ## clean the RAM
  gc()
  ##check for findings
  if(nrow(query_block_raw_df) == 0) {
    return(0)
  } else {
    ## build prot to genebank hash
    prot_2_genebank_db <- select(par_list["prot_info"], c(prot_id, genebank_id))
    prot_2_genebank_sub <- collect(filter(prot_2_genebank_db, prot_id %in% query_block_raw_df$prot_id))
    prot_2_genebank_map <- hashmap(prot_2_genebank_sub$prot_id, prot_2_genebank_sub$genebank_id)
    query_block_raw_df$genebank_id <- prot_2_genebank_map[[query_block_raw_df$prot_id]] 
    query_block_df <- filter(query_block_raw_df, pos_start < pos_end)
    return(query_block_df)
  }
}


