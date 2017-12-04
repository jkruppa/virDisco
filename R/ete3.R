##' Test
##'
##' Test
##' @title Test 
##' @param tax_ids 
##' @param out_file 
##' @param pdf_file 
##' @param clean 
##' @return NULL
##' @author Jochen Kruppa
##' @export
ete2_plot <- function(tax_ids, out_file, pdf_file, clean = FALSE){
  ## check for file
  if(nchar(file_ext(out_file)) != 0) stop("out_file should have no extension like '.tax'")
  tax_ids_full_file <- str_c(out_file, "_full.tax")
  write.table0(unique(tax_ids), tax_ids_full_file)
  ## check if all tax ids are available and drop the unkown
  tree_info_file <- str_c(out_file, "_info.tax")
  ete3_info_cmd <- paste("cut -f1", tax_ids_full_file,
                         "| ete3 ncbiquery --info",
                         ">", tree_info_file)
  runCMD(ete3_info_cmd)
  ## get the knwon tax ids into a tree 
  tree_info <- read.delim(tree_info_file)[,1]
  if(length(tree_info) > 1) {
    tax_ids_file <- str_c(out_file, "_sub.tax")
    write.table0(tree_info, tax_ids_file)
    ## run the plotting
    ete3_cmd <- paste("cut -f1", tax_ids_file,
                      "| ete3 ncbiquery --tree | xvfb-run ete3 view --ncbi --image",
                      pdf_file)
    runCMD(ete3_cmd)
  } else {
    talk("Only one species found. Tree makes no sense -> skip")
  }
  ## clean everything
  if(clean) unlink(list.files(dirname(out_file), ".tax$", full.names = TRUE))
}
## ------------------------------------------------------------
