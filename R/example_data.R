#' Covariance Model for 5.8S rRNA
#'
#' This is a convenience function to return the path the RFAM 5.8S rRNA
#' covariance model,
#' for use in [cmsearch()] or [cmalign()].
#' The original file is from <https://rfam.xfam.org/family/RF00002/cm>.
#'
#' @return (`character` string) the path to the CM file
#' @examples
#'     # search sequences from a fasta file for Rfam RF00002 (5.8S rRNA)
#'     cm <- cm_5_8S()
#'     sampfasta <- sample_rRNA_fasta()
#'     cmsearch(cm = cm, seq = sampfasta, cpu = 1)
#' @export
cm_5_8S <- function() {
    system.file(file.path("extdata", "RF00002.cm"), package = "inferrnal")
}

#' Annotated Seed Alignment for 5.8S rRNA
#'
#' This is a convenience function to return the path the seed alignment for the
#' RFAM 5.8S rRNA covariance model,
#' for use in [cmbuild()].
#' The original file is from
#' <https://rfam.xfam.org/family/RF00002/stockholm>.
#'
#' @return (`character` string) the path to the stk file
#' @export
#' @examples
#'     cmbuild(
#'         msafile = stk_5_8S(),
#'         cmfile_out = "/dev/null",
#'         force = TRUE,
#'         quiet = FALSE
#'     )
stk_5_8S <- function() {
    system.file(file.path("extdata", "RF00002.stk"), package = "inferrnal")
}

#' Sample rRNA sequences
#'
#' These are example files based on amplicon sequences from a soil metabarcoding
#' study. `sample_rRNA_fasta` is a FASTA file of amplicon sequences,
#' `sample_rRNA_5_8S` a FASTA file of the same amplicons, but including
#' only the 5.8S region, and `sample_rRNA_stk` is a Stockholm-format
#' alignment of the 5.8S region.
#'
#' @rdname sample_data
#' @return a path to the sample data file
#' @export
#' @examples
#'     # search sequences from a fasta file for Rfam RF00002 (5.8S rRNA)
#'     cm <- cm_5_8S()
#'     sampfasta <- sample_rRNA_fasta()
#'     cmsearch(cm = cm, seq = sampfasta, cpu = 1)
#'
#'     # read a stockholm alignment
#'     msafile <- sample_rRNA_stk()
#'     msa <- read_stockholm_msa(msafile)
#'     msa
sample_rRNA_fasta <- function() {
    system.file(file.path("extdata", "sample.fasta"), package = "inferrnal")
}

#' @rdname sample_data
#' @export
sample_rRNA_stk <- function() {
    system.file(file.path("extdata", "sample.stk"), package = "inferrnal")
}

#' @rdname sample_data
#' @export
sample_rRNA_5_8S <- function() {
    system.file(file.path("extdata", "samp_5_8S.fasta"), package = "inferrnal")
}
