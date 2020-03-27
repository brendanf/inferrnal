#' Covariance Model for 5.8S rRNA
#'
#' This is a convenience function to return the path the RFAM 5.8S rRNA covariance model,
#' for use in \code{\link{cmsearch}} or \code{\link{cmalign}}.
#' The original file is from \href{https://rfam.xfam.org/family/RF00002/cm}.
#'
#' @value (\code{character} string) the path to the CM file
#' @export
cm_5_8S <- function() system.file(file.path("extdata", "RF00002.cm"), package = "inferrnal")

#' Sample rRNA sequences
#'
#' This is a FASTA file of amplicon sequences
#' @export
sample_rRNA_fasta <- function() system.file(file.path("extdata", "sample.fasta"), package = "inferrnal")

#' @export
sample_rRNA_stk <- function() system.file(file.path("extdata", "sample.stk"), package = "inferrnal")

#' @export
sample_rRNA_5_8S <- function() system.file(file.path("extdata", "samp_5_8S.fasta"), package = "inferrnal")
