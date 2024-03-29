% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/example_data.R
\name{sample_rRNA_fasta}
\alias{sample_rRNA_fasta}
\alias{sample_rRNA_stk}
\alias{sample_rRNA_5_8S}
\title{Sample rRNA sequences}
\usage{
sample_rRNA_fasta()

sample_rRNA_stk()

sample_rRNA_5_8S()
}
\value{
a path to the sample data file
}
\description{
These are example files based on amplicon sequences from a soil metabarcoding
study. \code{sample_rRNA_fasta} is a FASTA file of amplicon sequences,
\code{sample_rRNA_5_8S} a FASTA file of the same amplicons, but including
only the 5.8S region, and \code{sample_rRNA_stk} is a Stockholm-format
alignment of the 5.8S region.
}
\examples{
    # search sequences from a fasta file for Rfam RF00002 (5.8S rRNA)
    cm <- cm_5_8S()
    sampfasta <- sample_rRNA_fasta()
    cmsearch(cm = cm, seq = sampfasta, cpu = 1)

    # read a stockholm alignment
    msafile <- sample_rRNA_stk()
    msa <- read_stockholm_msa(msafile)
    msa
}
