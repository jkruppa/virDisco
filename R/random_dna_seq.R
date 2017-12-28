##' descr
##'
##' details
##' @title title
##' @param inFile 
##' @param outFile 
##' @param k 
##' @param n 
##' @param clean 
##' @param tmpDir 
##' @param fasta_formatter 
##' @param fasta_uShuffle 
##' @return Null
##' @author Jochen Kruppa
##' @export
random_dna_seq <- function(inFile, outFile, k, n = 1, clean = TRUE,
                           tmpDir, fasta_formatter = NULL,
                           fasta_uShuffle = NULL){
  ## both needed programs
  if(is.null(fasta_formatter)){
    fasta_formatter <- file.path(programDir, "fastx-Toolkit/bin/fasta_formatter")
  }
  if(is.null(fasta_uShuffle)) {
    fasta_uShuffle <- file.path(programDir, "fasta_ushuffle/fasta_ushuffle")
  }
  ## build up tmpFile
  tmpFile <- file.path(tmpDir, gsub(".fasta|.fa", "_tmp.shuffle", basename(inFile)))
  ## build up commands: first fasta with no line breaks
  formatterCMD <- paste(fasta_formatter,
                        "-i", inFile,
                        "-o", tmpFile)
  try(system(formatterCMD))
  ## build up commands: second shuffle the tmp file
  uShuffleCMD <- paste(fasta_uShuffle,
                       "<", tmpFile,
                       "-n", n,
                       "-k", k,
                       "-r", 1000,
                       ## "-seed", 20160916,
                       ">", outFile,
                       "2>", file.path(tmpDir, "shuffle.log"))
  try(system(uShuffleCMD))
  ## clean...
  if(clean){
    file.remove(tmpFile)
  }
}
