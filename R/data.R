#' Output of ord_findings
#'
#' A example output of \code{\link{ord_findings}}.
#'
#' @format A data.frame with five entries:
#' \describe{
#'   \item{genebank_id}{Genbank ID}
#'   \item{dna_reads}{Count of DNA reads mapped to the reference Genbank ID}
#'   \item{dna_rank}{Rank of the DNA read count}
#'   \item{amino_reads}{Count of AA reads mapped to the reference Genbank ID}
#'   \item{amino_rank}{Rank of the AA read count}
#' }
#' @source \url{ftp://ftp.ncbi.nlm.nih.gov/genbank}
"NHS_10001_count_tbl"

#' Output of map_pep_ref
#'
#' A example output of \code{\link{map_pep_ref}}.
#'
#' @format A list with two data.frame:
#' $count_df
#' \describe{
#'   \item{genebank_id}{Genbank reference ID}
#'   \item{count_prot}{Number of translated reads mapped to the genebank ID reference} 
#' }
#' $alignment_df
#' \describe{
#'   \item{qname}{Read name i.e. query name}
#'   \item{prot_id}{Protein reference ID}
#'   \item{pos_start}{Poistion start on the reference DNA sequence}
#'   \item{pos_end}{Poistion end on the reference DNA sequence}
#'   \item{score_bit}{Mapping quality}
#'   \item{score_expected}{e-value of Blast}
#'   \item{ident}{How equal is the reference to the read?}
#'   \item{genebank_id}{Genebank ID of the reference sequence}
#' }
#' @source \url{ftp://ftp.ncbi.nlm.nih.gov/genbank}
"NHS_10001_map_aa_list"


#' Output of get_consensus_df
#'
#' A example output of the \code{\link{get_consensus_df}} function.
#'
#' @format A list with two entries:
#'
#' $consensus_df
#' \describe{
#'   \item{genebank_id}{Genbank ID of the reference viral strain}
#'   \item{con_prob}{Coverage probability. How unique is the base?}
#'   \item{pos}{Position on the reference}
#'   \item{bool}{Mathc between the reference and the consensus sequence}
#' }
#'
#' $entropy_df
#' \describe{
#'   \item{genebank_id}{Genbank ID of the reference viral strain}
#'   \item{entropy}{Entropy of the consensus reference}
#' }
#'
#' $run_len_df
#' \describe{
#'   \item{genebank_id}{Genbank ID of the reference viral strain}
#'   \item{min}{Minimum run of bases to the reference}
#'   \item{first}{First quartile of the runs to the reference}
#'   \item{median}{Median of the runs to the reference}
#'   \item{mean}{Mean of the runs to the reference}
#'   \item{third}{Third quartile of the runs to the reference}
#'   \item{max}{Maximum run of bases to the reference}
#'   \item{run_to_ref}{Maximal run to the reference}
#'   \item{run_not_to_ref}{Maximal run _not_ to the reference}
#' }
#'
#' $median_con_prob_df
#' \describe{
#'   \item{genebank_id}{Genbank ID of the reference viral strain}
#'   \item{min}{Minimum of the consensus_df$con_prob}
#'   \item{first}{First quartile of the consensus_df$con_prob}
#'   \item{median}{Median of the consensus_df$con_prob}
#'   \item{mean}{Mean of the consensus_df$con_prob}
#'   \item{third}{Third quartile of the consensus_df$con_prob}
#'   \item{max}{Maximum of the consensus_df$con_prob}
#' }
#'
#'$coverage_equal_ref_df
#' \describe{
#'   \item{genebank_id}{Genbank ID of the reference viral strain}
#'   \item{coverage_equal_ref}{How many of the consensus base match the reference base?}
#' }
#'
#' $coverage_overall_ref_df
#' \describe{
#'   \item{genebank_id}{Genbank ID of the reference viral strain}
#'   \item{coverage_overall_ref}{Ho much of the reference is covered by the consensus sequence}
#' }
#' @source \url{ftp://ftp.ncbi.nlm.nih.gov/genbank}
"NHS_10001_consensus"



#' Output of map_dna_ref
#'
#' A example output of \code{\link{map_dna_ref}}.
#'
#' @format A list with two data.frame:
#' 
#' $count_df
#' \describe{
#'   \item{genebank_id}{Genbank reference ID}
#'   \item{count_read}{Number of reads mapped to the genebank ID reference}
#' }
#' 
#' $alignment_df
#' \describe{
#'   \item{genebank_id}{Genbank reference ID}
#'   \item{qname}{Read name i.e. query name}
#'   \item{strand}{Strand}
#'   \item{read_length}{Read length}
#'   \item{pos_start}{Poistion start on the reference DNA sequence}
#'   \item{pos_end}{Poistion end on the reference DNA sequence}
#'   \item{mapq}{Mapping quality}
#'   \item{seq_id}{Additional identifier}
#' }
#' @source \url{ftp://ftp.ncbi.nlm.nih.gov/genbank}
"NHS_10001_map_dna_list"


#' Example fastq files
#'
#' A list with two fastq objects including R1 and R2 reads
#'
#' @format A list with two entries:
#' \describe{
#'   \item{R1}{Fastq reads object with 10000 paired reads of R1}
#'   \item{R2}{Fastq reads object with 10000 paired reads of R2}
#' }
#' @source \url{ftp://ftp.ncbi.nlm.nih.gov/genbank}
"infected_fastq"


#' Example species infro db
#'
#' A dataset containing 287 viral strains and information
#'
#' @format A dataset containing 287 viral strains and information:
#' \describe{
#'   \item{1st}{Genbank ID}
#'   \item{2nd}{Molecular type} 
#'   \item{3rd}{Strain: RNA, DNA or other}
#'   \item{4th}{Taxonomic ID}
#'   \item{5th}{Additional accession number}
#'   \item{6th}{Length}
#'   \item{7th}{Organism}
#'   \item{8th}{Description}
#' }
#' @source \url{ftp://ftp.ncbi.nlm.nih.gov/genbank}
"strain_info_db"

#' Example amino acid information db
#'
#' A dataset containing 519 protein IDs, the start and the end position
#' on the DNA sequence
#'
#' @format A data.frame with 519 protein IDs and the connected genbank IDs. 
#' \describe{
#'   \item{1st}{Protein ID}
#'   \item{2nd}{Genbank ID}
#'   \item{3rd}{Protein start position}
#'   \item{4th}{Protein end position}
#'   \item{5th}{Identifier}
#' }
#' @source \url{ftp://ftp.ncbi.nlm.nih.gov/genbank}
"aa_info_db"

#' DNAStringSet of DNA sequences
#'
#' A DNAStringSet containing 287 DNA sequences from the NCBI GenBank.
#'
#' @format A named DNAStringSet with 287 viral entries. The names are split by ' ' meaning:
#' \describe{
#'   \item{1st}{Genbank ID}
#'   \item{2nd}{Genbank ID}
#'   \item{3rd - End}{Description}
#' }
#' @source \url{ftp://ftp.ncbi.nlm.nih.gov/genbank}
"dna_seqs"


#' AAStringSet of amino acid sequences
#'
#' A dataset containing 521 amino acid sequences from the
#' corresponding DNA sequencing fasta file
#'
#' @format A AAStringSet with 521 viral gene entries. The names are split by '_' meaning:
#' \describe{
#'   \item{1st}{Protein ID}
#'   \item{2nd}{Genbank ID}
#'   \item{3rd}{Starting base on DNA sequence}
#'   \item{4th}{Ending base on DNA sequence}
#' }
#' @source \url{ftp://ftp.ncbi.nlm.nih.gov/genbank}
"aa_seqs"


