##' Test
##'
##' Test
##' @title Test 
##' @param inFile 
##' @param samOutFile 
##' @param referenceDir 
##' @param method 
##' @param p 
##' @param all 
##' @return character
##' @author Jochen Kruppa
##' @export
get_bowtie2_cmd <- function(inFile, samOutFile, referenceDir, method, p = 32, all = FALSE){
  require(tools)
  ## get the params for bowtie together
  params <- vector(mode = "character")
  ## check for fasta
  if(any(file_ext(inFile) %in% c("fasta", "fa"))){
    params <- paste("-f")
  }
  ## check for paired or single
  if(length(unlist(inFile)) == 2){
    params <- paste(params,
                    "-1", unlist(inFile)[1],
                    "-2", unlist(inFile)[2],
                    "--no-mixed")
  } else {
    params <- paste(params, "-U", inFile)
  }   
  params <- paste(params,
                  switch(method,
                         "end-to-end" = {
                           paste(
                             "--very-fast"## ,
                             ## str_c("--score-min L,0,-", 0.95)
                           )
                         },
                         "local" = {
                           paste(
                             "--very-fast-local",                                 
                             str_c("--score-min G,1,", 1))
                         })
                  )
  ## get all reads
  if(!all) {
    params <- paste(params, "--no-unal")
  }
  ## build bowtie2 command with params
  bowtie2CMD <- paste("bowtie2",
                      paste("-x", referenceDir,
                            "-p", p, ## number of cores
                            "-k", 10, ## allow 10 hits per read
                            ## "-a", ## find all aligments
                            ## "--no-unal", ## remove all unaligned reads from sam
                            params, 
                            "-S", samOutFile,
                            "--quiet")
                      )
  return(bowtie2CMD)
}
