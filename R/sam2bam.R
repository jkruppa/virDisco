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
    talk("[SAMTOOLS] Working on ", basename(x))
    if(index_sort){
    runCMD(paste(program_list["samtools"], "view -b -S", x, ">",
                 gsub(".sam", "_tmp.bam", x),
                 "2>", file.path(tmpDir, "samtools.log"))
           )
    runCMD(paste(program_list["samtools"], "sort",
                 gsub(".sam", "_tmp.bam", x), ">", gsub(".sam", ".bam", x),
                 "2>", file.path(tmpDir, "samtools.log"))
           )
    runCMD(paste(program_list["samtools"], "index", gsub(".sam", ".bam", x),
                 "2>", file.path(tmpDir, "samtools.log"))
           )
    } else {
      runCMD(paste(program_list["samtools"], "view -b -S", x, ">",
                   gsub(".sam", ".bam", x),
                   "2>", file.path(tmpDir, "samtools.log"))
             )
    }
    if(!keep){
      unlink(c(x,
               gsub(".sam", "_tmp.bam", x),
               gsub(".sam", "_tmp.bam", x)
               ))
    }
  })
}
