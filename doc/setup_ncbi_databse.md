# Setup the NCBI GenBank database

In the following the download and processing of the NCBI GenBank *viral* database is demonstrated. We use hier chunks of R code and the user might like to adjust the download.

### File setup for the NCBI GenBank database

While the run of the virDisco many files are produced and needed. Therefore a good file system managment should be used and setup. It is not a good idea to store and save all the input and output files into one single folder. 

All the files will be stored in the `main_dir`, which is located in the tmp folder of your PC. Change this line to use a different folder as root folder. Further, we need a tmp folder, here `tmp_viral_dir`,  to store all the temporal files. In this example it does not really matter but in a real example the intermediate files become really big and must not be stored.

```R
library(pacman) ## install pacman if needed
p_load(RCurl, stringr, plyr)

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

Normally over 50 files will be downloaded and further prcessed. To speedup everything, we only select four files.

```R
## CAUTION we only select four files!
genbank_file_virus <- genbank_file_virus[1:4] 
```

All downloads are done only from now on with a fraction of all GenBank databasse files!

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
