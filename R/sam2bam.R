##' A serial function to convert sam files from the mapping into bam files.
##'
##' This is not a standalone function. The samtools path is needed in
##' the program_list(). 
##' @title Wrapper function for samtools
##' @param sam_files A vector of sam_files to be converted
##' @param keep Should the sam files be keeped? [default = FALSE]
##' @param index_sort Should the bam file sbe sorted? [default =
##'   FALSE]
##' @return NULL
##' @author Jochen Kruppa
##' @export
sam2bam_serial <- function(sam_files, keep = FALSE, index_sort = FALSE){
  l_ply(sam_files, function(x){
    talk("[SAMTOOLS] Working on ", basename(x))
    if(index_sort){
    runCMD(paste(program_list["samtools"], "view -b -S", x, ">",
                 gsub(".sam", "_tmp.bam", x),
                 "2>", file.path(par_list["tmp_dir"], "samtools.log"))
           )
    runCMD(paste(program_list["samtools"], "sort",
                 gsub(".sam", "_tmp.bam", x), ">", gsub(".sam", ".bam", x),
                 "2>", file.path(par_list["tmp_dir"], "samtools.log"))
           )
    runCMD(paste(program_list["samtools"], "index", gsub(".sam", ".bam", x),
                 "2>", file.path(par_list["tmp_dir"], "samtools.log"))
           )
    } else {
      runCMD(paste(program_list["samtools"], "view -b -S", x, ">",
                   gsub(".sam", ".bam", x),
                   "2>", file.path(par_list["tmp_dir"], "samtools.log"))
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
