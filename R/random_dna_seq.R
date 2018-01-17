##' The function shuffles a given sequence and keeps a given k
##' distribution
##'
##' The function calls uShuffle, which can only handle fasta files
##' without a linebreak. Therefore, the program formatter removes all
##' linebreaks from the fasta file beforehand. Than uShuffle shuffles
##' the sequence by keeping a distrubtion of kmer's given k.
##' @title Shuffle a sequence with a kept kmer distribution
##' @param inFile File path to the fasta file to be shuffled
##' @param outFile File path to the out dir
##' @param k Kept k for the shuffling (k < 25)
##' @param n Number of produced shuffle sequences per sequence
##'   [default = 1]
##' @param clean Should everthing removed [default = TRUE]
##' @param tmpDir Where to save temp files
##' @param fasta_formatter File path to the fasta formatter
##'   executables
##' @param fasta_uShuffle File path to the fasta uShuffle executables
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
