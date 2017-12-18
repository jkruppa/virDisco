##' Test
##'
##' Test
##' @title Test 
##' @param sam_files 
##' @param keep 
##' @param index_sort 
##' @return NULL
##' @author Jochen Kruppa
##' @export
sam2bam_serial <- function(sam_files, keep = FALSE, index_sort = FALSE){
  l_ply(sam_files, function(x){
    talk("Working on ", basename(x))
    if(index_sort){
    runCMD(paste(program_list["samtools"], "view -b -S", x, ">",
                 gsub(".sam", "_tmp.bam", x)))
    runCMD(paste(program_list["samtools"], "sort",
                 gsub(".sam", "_tmp.bam", x), ">", gsub(".sam", ".bam", x)))
    runCMD(paste(program_list["samtools"], "index", gsub(".sam", ".bam", x)))
    } else {
      runCMD(paste(program_list["samtools"], "view -b -S", x, ">",
                   gsub(".sam", ".bam", x)))      
    }
    if(!keep){
      unlink(c(x,
               gsub(".sam", "_tmp.bam", x),
               gsub(".sam", "_tmp.bam", x)
               ))
    }
  })
}
