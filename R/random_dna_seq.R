##' descr
##'
##' details
##' @title title
##' @param inFile 
##' @param outFile 
##' @param k 
##' @param n 
##' @param clean 
##' @return Null
##' @author Jochen Kruppa
##' @export
random_dna_seq <- function(inFile, outFile, k, n = 1, clean = TRUE){
    ## both needed programs
    fasta_formatter <- file.path(programDir, "fastx-Toolkit/bin/fasta_formatter")
    fasta_uShuffle <- file.path(programDir, "fasta_ushuffle/fasta_ushuffle")
    ## build up tmpFile
    tmpFile <- file.path(tmpDir, gsub(".fasta|.fa", "_tmp.shuffle", basename(inFile)))
    ## build up commands: first fasta with no line breaks
    formatterCMD <- paste(fasta_formatter,
                          "-i", inFile,
                          "-o", tmpFile)
    runCMD(formatterCMD)
    ## build up commands: second shuffle the tmp file
    uShuffleCMD <- paste(fasta_uShuffle,
                         "<", tmpFile,
                         "-n", n,
                         "-k", k,
                         "-r", 1000,
                         ## "-seed", 20160916,
                         ">", outFile,
                         "2>", file.path(tmpDir, "shuffle.log"))
    runCMD(uShuffleCMD)
    ## clean...
    if(clean){
        file.remove(tmpFile)
    }
}
