#' Covariance Model for 5.8S rRNA
#'
#' This is a convenience function to return the path the RFAM 5.8S rRNA covariance model,
#' for use in \code{\link{cmsearch}} or \code{\link{cmalign}}.
#' The original file is from \url{https://rfam.xfam.org/family/RF00002/cm}.
#'
#' @return (\code{character} string) the path to the CM file
#' @export
cm_5_8S <- function() {
    system.file(file.path("extdata", "RF00002.cm"), package = "inferrnal")
}

#' Annotated Seed Alignment for 5.8S rRNA
#'
#' This is a convenience function to return the path the seed alignment for the
#' RFAM 5.8S rRNA covariance model,
#' for use in \code{\link{cmbuild}}.
#' The original file is from
#' \url{https://rfam.xfam.org/family/RF00002/stockholm}.
#'
#' @return (\code{character} string) the path to the stk file
#' @export
stk_5_8S <- function() {
    system.file(file.path("extdata", "RF00002.stk"), package = "inferrnal")
}

#' Sample rRNA sequences
#'
#' These are example files based on amplicon sequences from a soil metabarcoding
#' study. \code{sample_rRNA_fasta} is a FASTA file of amplicon sequences,
#' \code{sample_rRNA_5_8S} a FASTA file of the same amplicons, but including
#' only the 5.8S region, and \code{sample_rRNA_stk} is a Stockholm-format
#' alignment of the 5.8S region.
#'
#' @rdname sample_data
#' @return a path to the sample data file
#' @export
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
