#' Dataset
#'
#' A dataset containing 1000 DNA sequences from the NCBI GenBank.
#'
#' @format A fasta file with 1000 entries. The names are split by ' ' meaning:
#' \describe{
#'   \item{1st}{Genbank ID}
#'   \item{2nd}{Genbank ID}
#'   \item{3rd - End}{Description}
#' }
#' @source \url{ftp://ftp.ncbi.nlm.nih.gov/genbank}
"strain_info_db"

#' Dataset
#'
#' A dataset containing 1000 DNA sequences from the NCBI GenBank.
#'
#' @format A fasta file with 1000 entries. The names are split by ' ' meaning:
#' \describe{
#'   \item{1st}{Genbank ID}
#'   \item{2nd}{Genbank ID}
#'   \item{3rd - End}{Description}
#' }
#' @source \url{ftp://ftp.ncbi.nlm.nih.gov/genbank}
"aa_info_db"

#' DNAStringSet of DNA sequences
#'
#' A dataset containing 1000 DNA sequences from the NCBI GenBank.
#'
#' @format A named DNAStringSet with 1000 viral entries. The names are split by ' ' meaning:
#' \describe{
#'   \item{1st}{Genbank ID}
#'   \item{2nd}{Genbank ID}
#'   \item{3rd - End}{Description}
#' }
#' @source \url{ftp://ftp.ncbi.nlm.nih.gov/genbank}
"dna_seqs"


#' AAStringSet of amino acid sequences
#'
#' A dataset containing 1268 amino acid sequences from the
#' corresponding DNA sequencing fasta file
#'
#' @format A AAStringSet with 1268 viral gene entries. The names are split by '_' meaning:
#' \describe{
#'   \item{1st}{Protein ID}
#'   \item{2nd}{Genbank ID}
#'   \item{3rd}{Starting base on DNA sequence}
#'   \item{4th}{Endinf base on DNA sequence}
#' }
#' @source \url{ftp://ftp.ncbi.nlm.nih.gov/genbank}
"aa_seqs"


