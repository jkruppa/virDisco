% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trimming.R
\name{fastq_trimmer}
\alias{fastq_trimmer}
\title{Internal trimmomatic call function}
\usage{
fastq_trimmer(inFile, outFile, leading = 10, trailing = 10,
  minlength = 50, illuminaclip, log_file)
}
\arguments{
\item{inFile}{Infile of the sample, see \code{\link{get_sub_files}}.}

\item{outFile}{Outfile of the sample, see \code{\link{get_sub_files}}.}

\item{leading}{See trimmomatic manual}

\item{trailing}{See trimmomatic manual}

\item{minlength}{See trimmomatic manual}

\item{illuminaclip}{See trimmomatic manual}

\item{log_file}{file path to log file}
}
\value{
Test
}
\description{
Internal trimmomatic call function
}
\details{
Internal trimmomatic call function
}
\author{
Jochen Kruppa
}
