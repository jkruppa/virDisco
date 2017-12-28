# Setup the NCBI GenBank database

In the following the download and processing of the NCBI GenBank *viral* database is demonstrated. We use hier chunks of R code and the user might like to adjust the download.

### File setup for the NCBI GenBank database

While the run of the virDisco many files are produced and needed. Therefore a good file system managment should be used and setup. It is not a good idea to store and save all the input and output files into one single folder. 

All the files will be stored in the `main_dir`, which is located in the tmp folder of your PC. Change this line to use a different folder as root folder. Further, we need a tmp folder, here `tmp_viral_dir`,  to store all the temporal files. In this example it does not really matter but in a real example the intermediate files become really big and must not be stored.

```R
library(pacman) ## install pacman if needed
p_load(RCurl, stringr, plyr, Biostrings)

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

In the first step, we will download all ftp GenBank files connected to a virus. This is indicated by `gbvrl.*` at the beginning of  a database file. We wait 5 seconds after each download to avoid any ftp timeout. 

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

Normally over 50 files will be downloaded and further processed. To speedup everything, we only select four files.

```R
## CAUTION we only select four files!
genbank_file_virus <- genbank_file_virus[1:4] 
```

All steps are done from now only with a fraction of all GenBank databasse files! Remove the above code line, to process all files.

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

## Step 2.1: Add the decoy database

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






