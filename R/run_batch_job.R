##' Control function for the processing of many fastq samples of one folder
##'
##' This is a main control function to get all the fastq files of one
##' folder processed. The output is written in a other folder
##' separated by the file names
##' @title Function to run the pipeline on one folder
##' @param batch_in_dir Path to the folder with the fastq files
##' @param batch_out_dir Path to the ouput folder
##' @param par_list See \code{link{set_par_list}} for the possible
##'   parameters
##' @param run_samples Character vector. Should special samples be run
##'   and not all? [default = NULL]
##' @return NULL
##' @author Jochen Kruppa
##' @export
run_batch_job <- function(batch_in_dir, batch_out_dir, par_list, run_samples = NULL) {
  ## get the batch in files
  batch_file_list <- prepare_file_list(batch_in_dir)
  ## get only some samples
  if(!is.null(run_samples)) {
    if(all(run_samples %in% names(batch_file_list))){
      batch_file_list <- batch_file_list[names(batch_file_list) %in% run_samples]
    } else {
      stop("Cannnot find the run_samples in names(batch_file_list)!")
    }
  }
  ## start the processing of the files
  l_ply(seq_along(batch_file_list), function(i){
    sample <- names(batch_file_list)[i]
    message("+------------------------------------------------------+")
    message("|                                                      |")
    message(str_c("| Start sample: ", sample, "                              |"))
    message("|                                                      |")
    message("+------------------------------------------------------+")
    sample_dir <- file.path(batch_out_dir, sample)
    if(!file.exists(sample_dir)) dir.create(sample_dir)
    sample_out <- file.path(sample_dir, sample)
    par_list["pdf_file"] <- file.path(sample_dir, str_c(sample, ".pdf"))
    artificial_genome_mapping(in_file = batch_file_list[i],
                              out_file = sample_out,
                              par_list
                              )
  })
  message("+------------------------------------------------------+")
  talk("Copy and zip the pdf and xlsx files")
  message("+------------------------------------------------------+")
  ## copy everthing back and pack it
  batch_pdf_files <- grep("contig", dir0(batch_out_dir, pattern = ".pdf$"), invert = TRUE, value = TRUE)
  batch_xlsx_files <- dir0(batch_out_dir, pattern = ".xlsx$")
  if(all(grepl("TIHO", batch_file_list))) {
    batch_res_name <- unique(gsub(".*(TIHO.*)\\/.*", "\\1", batch_file_list))
  } else {
    batch_res_name <- basename(dirname(batch_in_dir))
  }
  batch_res_dir <- file.path(dirname(batch_out_dir),
                             str_c(batch_res_name, "_results"))
  message(batch_res_dir)
  if(!file.exists(batch_res_dir)){
    dir.create(batch_res_dir)
  } else {
    batch_res_dir <- str_c(batch_res_dir, "_", stringi::stri_rand_strings(1, 10)) 
    dir.create(batch_res_dir)
  }
  ## copy pdf and xlsx 
  file.copy(batch_pdf_files, batch_res_dir)
  file.copy(batch_xlsx_files, batch_res_dir)
  ## zip everything
  runCMD(paste("zip -r -j -q", batch_res_dir, batch_res_dir))
  message("+------------------------------------------------------+")
  talk("Finished")
  message("+------------------------------------------------------+") 
}
