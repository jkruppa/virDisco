# Setup the NCBI GenBank database

In the following the download and processing of the NCBI GenBank *viral* database is demonstrated. We use hier chunks of R code and the user might like to adjust the download.

*Please read first the full tutorial. Some programs must be installed and it might be feasiable to adjust some code chunks for your own purpose. Overall 2.4 million sequences will be parsed. Therefore, start with a small amount of files from NCBI and see what you get.*

## Table of Contents
1. [File setup for the NCBI GenBank database](#file-setup-for-the-ncbi-genbank-database)
2. [Step 1: Download all viral database files](#step-1-download-all-viral-database-files)
3. [Step 2: Extract all sequences and feature information](#step-2-extract-all-sequences-and-feature-information)
4. [Step 3: Add the decoy database](#step-3-add-the-decoy-database)
5. [Step 4.1: Build bowtie-index on DNA data](#step-41-build-bowtie-index-on-dna-data)
6. [Step 4.2: Build star-index on DNA data](#step-42-build-star-index-on-dna-data)
7. [Step 5: Extract all peptide sequences](#step-5-extract-all-peptide-sequences)
8. [Step 6: Build pauda-index](#step-6-build-pauda-index)
9. [Step 7.1: SQlite database of gene information](#step-71-sqlite-database-of-gene-information)
10. [Step 7.2: SQlite database of species and description information](#step-72-sqlite-database-of-species-and-description-information)

## File setup for the NCBI GenBank database

While the run of the virDisco many files are produced and needed. Therefore a good file system managment should be used and setup. It is not a good idea to store and save all the input and output files into one single folder. 

Here, all the files will be stored in the `main_dir`, which is located in the tmp folder of your PC. Change this line to use a different folder as root folder. Further, we need a tmp folder, here `tmp_viral_dir`,  to store all the temporal files. In this example it does not really matter but in a real example the intermediate files become really big and must not be stored.

```R
library(pacman) ## install pacman if needed
p_load(RCurl, stringr, plyr, Biostrings, readr, tidyverse, ShortRead)

## setup the file system in the temp dir
main_dir <- tempdir()
tmp_viral_dir <- file.path(main_dir, "tmp")
sql_viral_dir <- file.path(main_dir, "sql")
genbank_ncbi_dir <- file.path(main_dir, "genbank_ncbi", "viral")
dir.create(genbank_ncbi_dir, recursive = TRUE)
dir.create(tmp_viral_dir, recursive = TRUE)
dir.create(sql_viral_dir, recursive = TRUE)

## Pauda dir, set to your own installation
paudaDir <- file.path("/home/programs/pauda-1.0.1/pauda-1.0.1/bin")
pauda_run <- file.path(paudaDir, "pauda-run")
pauda_build <- file.path(paudaDir, "pauda-build")
```

## Step 1: Download all viral database files

In the first step, we will download all NCBI GenBank files connected to a virus from the GenBank ftp server. The files are indicated by `gbvrl.*` for a virus at the beginning of each database file. We wait 5 seconds after each download to avoid any ftp timeout. Visit ftp://ftp.ncbi.nlm.nih.gov/genbank for a full overview of available data.

```R
## get all genbank files
genbank_file_names_raw <- getURL("ftp://ftp.ncbi.nlm.nih.gov/genbank/",
                                 verbose = TRUE, ftp.use.epsv = TRUE,
                                 dirlistonly = TRUE)

## process and select only viral strains
genbank_file_names <- str_split(genbank_file_names_raw, "\n", simplify = TRUE)[1,]
genbank_file_virus_names <- grep("gbvrl.*", genbank_file_names, value = TRUE)
genbank_file_virus <- str_c("ftp://ftp.ncbi.nlm.nih.gov/genbank/", genbank_file_virus_names)
```
Normally over *50 files* will be downloaded and further processed. To speedup everything, we only select *four files*.

```R
## CAUTION we only select four files!
genbank_file_virus <- genbank_file_virus[1:4] 
```

*All steps are done from now only with a fraction of all GenBank databasse files! Remove the above code line, to process all files.*

```R
## download all files to the genbank_ncbi_dir
setwd(genbank_ncbi_dir)
l_ply(genbank_file_virus, function(x) {
  if(any(str_detect(basename(x), dir(genbank_ncbi_dir)))){
    message("File ", basename(x), " already downloaded...")
  } else {
    message("Download ", basename(x))
    system(paste("wget -q", x))
    system(paste("gunzip", basename(x)))
    message("\tSleep 5 sec to avoid ftp timeout...")
    Sys.sleep(5) ## avoid ftp timeout
  }
})

## check if everything is there
all(str_replace(basename(genbank_file_virus), ".gz", "") %in% dir(genbank_ncbi_dir))

## get all the sequence files
genbank_ncbi_dna_files <- dir(genbank_ncbi_dir, pattern = ".seq$", full.names = TRUE)
```

Now, all files are available in the `embl` format. This is nice, because we have the sequence and the gene information. However, the format is not easy to parse. We use therefor the EMBOSS Tools (https://www.ebi.ac.uk/Tools/emboss/).

## Step 2: Extract all sequences and feature information

For the extraction of the sequence and protein information we use the EMBOSS Tools (https://www.ebi.ac.uk/Tools/emboss/). The function `seqret` and `extractfeat` allow to extract the sequence and amino information out of the seq files.  

```R
l_ply(genbank_ncbi_dna_files, function(x) {
  tmp_dna_file <- file.path(tmp_viral_dir, str_c(basename(x), "_dna.fa"))
  seqretCMD <- paste("seqret", 
                     "-sequence", x,
                     "-outseq", tmp_dna_file,
                     "-auto") # be quiet
  try(system(seqretCMD, wait = TRUE))
  tmp_cds_amino_file <- file.path(tmp_viral_dir, str_c(basename(x), "_cds_amino.fa"))
  extractfeatCMD <- paste("extractfeat", 
                          "-sequence", x,
                          "-outseq", tmp_cds_amino_file,
                          "-type", "CDS",
                          "-describe", "'db_xref | protein_id | translation'",
                          "-join",
                          "-auto",
                          "1>", file.path(tmp_viral_dir, "emboss.log"), "2>&1") # be quiet
  try(system(extractfeatCMD, wait = TRUE))
}, .progress = "text")
```

Finally, we must combine all the single fasta files onto one fasta. This can be done very quick by the `cat` command.

```R
## combine all fasta files
cat_cmd <- paste("find", tmp_viral_dir, "-name '*_dna.fa' -exec cat {} \\; >",
                 file.path(genbank_ncbi_dir, "gbvrl_multi.fa"))
try(system(cat_cmd, wait = TRUE))
```

## Step 3: Add the decoy database

We want to shuffle the original sequences of the viral database by a preserved `k`. Therefore, we use the program `uShuffle` (http://digital.cs.usu.edu/~mjiang/ushuffle/). However, `uShuffle` can only handle sequences of the maximal length of 9e6. Therefore, we build up chunks including seuquences of the summed up length smaller than 9e6. 

```R
gbvrl_multi_seq <- readDNAStringSet(file.path(genbank_ncbi_dir, "gbvrl_multi.fa"))

length_gbvrl_multi_seq <- sum(as.numeric(width(gbvrl_multi_seq)))

shuffle_threshold <- 9e6
gbvrl_multi_seq_list <- getSeqChunksByThreshold(gbvrl_multi_seq, shuffle_threshold)

## save the list
write_rds(gbvrl_multi_seq_list, file.path(genbank_ncbi_dir, "gbvrl_multi_seq_list.RDS"))

## read the list
gbvrl_multi_seq_list <- read_rds(file.path(genbank_ncbi_dir, "gbvrl_multi_seq_list.RDS"))
```

Finally, we must remove some non `ACGT` letters from the sequences. This will take some time.

```R
## combine the single sequences to one sequence
## remove letters, which are not ACGT
gbvrl_multi_comb_seq_list <- llply(gbvrl_multi_seq_list, function(x) {
  DNAStringSet(gsub("[^ACGT]", "", Reduce(c, x)))
}, .parallel = TRUE)
```

In the next step the sequences will be shuffled. Therefore, two programs are needed. First, the program `formatter` form the `fastx-Toolkit` (http://hannonlab.cshl.edu/fastx_toolkit/) to make the fasta files two liner, with no line break and the program `uShuffle` for the shuffling of the sequence.

```R
fasta_formatter <- file.path("/home/programs", "fastx-Toolkit/bin/fasta_formatter")
fasta_uShuffle <- file.path("/home/programs", "fasta_ushuffle/fasta_ushuffle")

k <- 15
l_ply(seq_along(gbvrl_multi_comb_seq_list), function(i){
  chunk_in <- file.path(tmp_viral_dir, str_c("chunk_in_", str_pad(i, 3, pad = "0"), ".fa"))
  chunk_out <- file.path(tmp_viral_dir, str_c("chunk_out_", str_pad(i, 3, pad = "0"), ".fa"))
  writeFasta(DNAStringSet(gbvrl_multi_comb_seq_list[[i]]), chunk_in)
  ## decoy random sequence 
  random_dna_seq(fasta_formatter = fasta_formatter,
                 fasta_uShuffle = fasta_uShuffle,
                 inFile = chunk_in,
                 outFile = chunk_out,
                 k = k,
                 tmpDir = tmp_viral_dir)
  unlink(chunk_in)
}, .progress = "text")
```
In the last step, we combine the shuffled chunk files and the original viral sequences into one final file `gbvrl_multi_decoy_15.fa`.

```R
## get everything together in one file
chunk_out_files <- dir(tmp_viral_dir, "chunk_out", full.names = TRUE)
decoy_k <- DNAStringSet(Reduce(c, llply(chunk_out_files, function(x){
  message(x)
  readDNAStringSet(x)}
  )))
names(decoy_k) <- rep("decoy_15", length(decoy_k))
decoy_k_file <- file.path(genbank_ncbi_dir, str_c("viral_decoy_k", k , ".fa"))
writeXStringSet(decoy_k, decoy_k_file)

## remove all chunk files
unlink(chunk_out_files)

## now we build up the reference
gbvrl_multi_decoy_k_seq <- c(gbvrl_multi_seq, decoy_k)
gbvrl_multi_decoy_k_seq_file <- file.path(genbank_ncbi_dir, "gbvrl_multi_decoy_15.fa")
writeXStringSet(gbvrl_multi_decoy_k_seq, gbvrl_multi_decoy_k_seq_file)
```

## Step 4.1: Build bowtie-index on DNA data

You can run the mapping in virDisco with the Bowtie2 mapper or the Star mapper. Both mapper need a indexed reference genome, which will be build in the following.

```R
referenceViralMultiFile <- file.path(genbank_ncbi_dir, "gbvrl_multi_decoy_15.fa")
referenceViralBowtieDir <- file.path(dirname(genbank_ncbi_dir),
                                     "viral_multi_bowtie/viral_multi_bowtie")
if(!file.exists(dirname(referenceViralBowtieDir))) dir.create(dirname(referenceViralBowtieDir))

bowtie2BuildCMD <- paste("bowtie2-build",
                         "-q",
                         referenceViralMultiFile,
                         referenceViralBowtieDir)
try(system(bowtie2BuildCMD))
```

## Step 4.2: Build star-index on DNA data

```R
referenceViralStarDir <- file.path(dirname(genbank_ncbi_dir), "viral_multi_star/viral_multi_star")

dir.create(dirname(referenceViralStarDir), showWarnings = FALSE)    
star_build_CMD <- paste(STAR,
                        "--runMode", "genomeGenerate",
                        "--genomeDir", dirname(referenceViralStarDir), 
                        "--genomeFastaFiles", referenceViralMultiFile,
                        "--genomeChrBinNbits", 10, 
                        "--runThreadN", par$nCores)
try(system(star_build_CMD))
```

## Step 5: Extract all peptide sequences



```R
l_ply(genbank_ncbi_dna_files, function(x) {
  sample <- basename(x)
  talk("Working on ", sample)
  tmp_cds_amino_file <- file.path(tmpDir, str_c(sample, "_cds_amino.fa"))
  tmp_cds_set <- readDNAStringSet(tmp_cds_amino_file)
  talk("Get the amino information and write to file")
  aa_file <- file.path(tmpDir, str_c(sample, "_aa.fa"))
  if(file.exists(aa_file)) unlink(aa_file)
  l_ply(names(tmp_cds_set), function(x) {
    if(grepl("translation", x)){
      tmp_aa_set <- AAStringSet(gsub(".*translation=\\\"(.*)\\\"\\).*", "\\1", x))
      ## get the sequence name
      genbank_id  <- str_split(x, " ", simplify = TRUE)[,1]
      prot_id <- gsub(".*protein_id=\\\"(.+?)\\\".*", "\\1", x)
      names(tmp_aa_set) <- str_c(prot_id, genbank_id, sep = "_")
      writeXStringSet(tmp_aa_set, aa_file, append = TRUE)
    }
  }, .progress = "text")
})
```

```R
## combine all fasta files
cat_cmd <- paste("find", tmpDir, "-name '*_aa.fa' -exec cat {} \\; >",
                 file.path(genbank_ncbi_dir, "gbvrl_multi_aa.fa"))
runCMD(cat_cmd)
```

## Step 6: Build pauda-index


```R
## build PAUDA index by bowtie2
pauda_build_dir <- file.path(tmp_viral_dir, "viral_multi_pauda_bowtie")

##                                        
runCMD(paste(pauda_build,
             file.path(genbank_ncbi_dir, "gbvrl_multi_aa.fa"),
             pauda_build_dir))

file.copy(list.files(pauda_build_dir, full.names = TRUE),
          file.path(dirname(genbank_ncbi_dir), "viral_multi_pauda_bowtie"),
          recursive = TRUE)
```

## Step 7.1: SQlite database of gene information

```R
ncbi_aa_fa <- file.path(genbank_ncbi_dir, "gbvrl_multi_aa.fa")
aa_info_db_file <- file.path(sql_viral_dir, "gbvrl_aa_info.sqlite3")


ncbi_aa_seq <- readAAStringSet(ncbi_aa_fa)

setup_aa_info_sqlite(ncbi_aa_seq, db_file = aa_info_db_file)
```

## Step 7.2: SQlite database of species and description information

```R
genbank_ncbi_info_df <- tbl_df(ldply(genbank_ncbi_dna_files, function(x) {
  sample <- basename(x)
  message("Working on ", sample)
  tmp_feat_file <- file.path(tmp_viral_dir, str_c(sample, "_feat.fa"))
  ##
  message("\tExtract feats")
  extractfeatCMD <- paste("extractfeat", 
                          "-sequence", x,
                          "-outseq", tmp_feat_file,
                          "-type", "source",
                          "-describe", "'db_xref | strain | mol_type'",
                          "-join",
                          "-auto",
                          "1>", file.path(tmp_viral_dir, "emboss.log"), "2>&1") # be quiet
  runCMD(extractfeatCMD)
  ## collect feature information
  tmp_feat_set <- readDNAStringSet(tmp_feat_file)
  tmp_aa_info <- ldply(names(tmp_feat_set), function(x) {
    info_raw <- str_split(x, " ", simplify = TRUE)
    genebank_id <- str_split(info_raw[,1], "_", simplify = TRUE)[1]
    mol_type <- gsub(".*mol_type=\\\"(.*)\\\", strain=.*", "\\1", x)
    strain <- gsub(".*strain=\\\"(.*)\\\", db.*", "\\1", x)
    tax_id <- gsub(".*taxon:(.*)\\\"\\).*", "\\1", x)
    return(data.frame(genebank_id, mol_type, strain, tax_id))
  })
  ## get sequence information
  tmp_info_file <- file.path(tmp_viral_dir, str_c(sample, "_info.txt"))
  ##
  message("\tExtract info")
  infoseqCMD <- paste("infoseq", 
                      "-sequence", x,
                      "-outfile", tmp_info_file,
                      "-nocolumns",
                      "-delimiter", "'\t'",
                      "-nodatabase",
                      "-nopgc",
                      "-notype",
                      "-nousa",
                      #"-noname",
                      "-auto",
                      "1>", file.path(tmp_viral_dir, "emboss.log"), "2>&1") # be quiet
  runCMD(infoseqCMD)
  tmp_seq_info <- tbl_df(read.delim(tmp_info_file, as.is = TRUE))
  tmp_info <- left_join(tbl_df(tmp_aa_info), tbl_df(tmp_seq_info),
                        by = c("genebank_id" = "Name"))
  message("Finished\n")
  return(tmp_info)
}))
```

```R
## copy everything to sqlite
genbank_ncbi_info_db_file <- file.path(sql_viral_dir, "gbvrl_info.sqlite3")
genbank_ncbi_info_db <- src_sqlite(genbank_ncbi_info_db_file, create = TRUE)
genbank_ncbi_info_sqlite <- copy_to(genbank_ncbi_info_db, genbank_ncbi_info_df, temporary = FALSE)
```R



