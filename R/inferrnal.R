
parse_stockholm_msa_chunk <- function(x, pos, acc) {

    gc <- stringr::str_match(x, "#=GC +([^ ]+) +(.+)")
    gc <- gc[,2:3, drop = FALSE]
    gc <- gc[stats::complete.cases(gc), , drop = FALSE]
    for (i in seq_len(nrow(gc))) {
        attr(acc, gc[i,1]) <- paste0(attr(acc, gc[i,1]), gc[i,2])
    }

    x <- stringr::str_match(x, "^(\\d+\\|)?([^#][^ ]*) +([^ ]+)$")
    x <- x[,3:4, drop = FALSE]
    x <- x[stats::complete.cases(x), , drop = FALSE]
    for (i in seq_len(nrow(x))) {
        if (x[i,1] %in% names(acc)) {
            acc[[x[i,1]]] <- paste0(acc[[x[i,1]]], x[i,2])
        } else {
            acc[[x[i,1]]] <- x[i,2]
        }
    }
    acc
}

#' Parse a Multiple Alignment in Stockholm Format
#'
#' Parses Stockholm-format multiple alignment files, including "\code{#=GC}"
#' lines.  Other annotations are ignored.
#'
#' @param stockholm (\code{character} scalar) Path to a file to parse
#'
#' @return a \code{list}.  The sequences themselves are element
#'     "\code{alignment}", and \code{#=GC} annotations are given as
#'     \code{character} strings; e.g., the "\code{#=GC SS_cons}" annotation is
#'     stored as an element named "\code{SS_cons}".
#' @export
#'
#' @examples
#'     msafile <- system.file(
#'         file.path("extdata", "sample.stk"),
#'         package = "inferrnal"
#'     )
#'     msa <- read_stockholm_msa(msafile)
#'     msa$alignment
#'     # consensus secondary structure
#'     msa$SS_cons
#'     # reference sequence
#'     msa$RF
read_stockholm_msa <- function(stockholm) {
    assertthat::assert_that((assertthat::is.string(stockholm) &&
                                file.exists(stockholm)) ||
                                methods::is(stockholm, "connection"))

    seqs <-
        readr::read_lines_chunked(
            stockholm,
            readr::AccumulateCallback$new(
                parse_stockholm_msa_chunk,
                acc = list()
            )
        )
    out <- attributes(seqs)
    out[["alignment"]] <- Biostrings::RNAMultipleAlignment(unlist(seqs))
    out
}

#' Search for Covariance Models (CM) in a set of sequences.
#'
#' This function calls "\code{cmsearch}" from
#' \href{http://eddylab.org/infernal/}{Infernal}.  Infernal must be installed.
#' Many parameters are not included (yet!), and the function is focused on
#' retrieving the hits table and, optionally, producing an alignment.
#'
#' @param cm (character of length 1) Path to the covariance model (.cm) file.
#'  The covariance model must include calibration data from running
#'  "\code{cmcalibrate}".
#' @param seq (filename, character vector,
#'   \code{\link[Biostrings]{XStringSet}}, or
#'   \code{\link[ShortRead]{ShortRead}}) Sequences to search with the CM.
#'   If a filename, the file can be of any type supported by Infernal.
#' @param glocal (logical of length 1) Whether to run the search in glocal mode
#'   (global with respect to the model, local with respect to the sequence).
#'   When \code{TRUE}, the search is faster, but will fail to find matches with
#'   only partially homologous structure.
#' @param alignment (filename) A file to save the aligned hits to.  If given,
#'   the alignment is saved in Stockholm format with annotations for secondary
#'   structure, posterior probablility, etc.
#' @param cpu (integer of length 1) The number of threads to use in the search.
#'
#' @return a \code{\link[tibble]{tibble}} with columns:
#'     \itemize{
#'         \item{target_name}{ (character) the name of the target sequence}
#'         \item{taget_accession}{(character) the target's accession number}
#'         \item{query_name}{(character) the name of the query CM}
#'         \item{query_accession}{(character) the query CM's accession number}
#'         \item{mdl}{(character) the model type ("cm" or "hmm")}
#'         \item{mdl_from}{(integer) the start location of the hit in the model}
#'         \item{mdl_to}{(integer) the end location of the hit in the model}
#'         \item{seq_from}{(integer) the start location of the hit in the
#'             sequence}
#'         \item{seq_to}{(integer) the end location of the hit in the sequence}
#'         \item{strand}{(character) the strand the hit was found on ("+" or
#'             "-")}
#'         \item{trunc)}{(character) whether the hit is truncated, and where
#'             ("no", "5'", "3'", "5'&3'", or "-" for hmm hits).}
#'         \item{pass}{(integer) which algorithm pass the hit was found on.}
#'         \item{gc}{(numeric) GC content of the hit}
#'         \item{bias}{(numeric) biased composition correction.  See the
#'             Infernal documentation.}
#'         \item{score}{(numeric) bit-score of the hit, including the biased
#'             composition correction.}
#'         \item{E_value}{(numeric) Expectation value for the hit.}}
#'         \item{inc}{(character) "!" if the sequence meets the inclusion
#'             threshold, "?" if it only meets the reporting threshold.}
#' @export
#'
#' @examples
#'     # search sequences from a fasta file for Rfam RF00002 (5.8S rRNA)
#'     cm5_8S <- system.file(
#'         file.path("extdata", "RF00002.cm"),
#'         package = "inferrnal"
#'     )
#'     sampfasta <- system.file(
#'         file.path("extdata", "sample.fasta"),
#'         package = "inferrnal"
#'     )
#'     cmsearch(cm = cm5_8S, seq = sampfasta, cpu = 1)
#'     # also works if the fasta file has already been loaded
#'     samp <- Biostrings::readDNAStringSet(sampfasta)
#'     cmsearch(cm = cm5_8S, seq = samp, cpu = 1)
cmsearch <- function(cm, seq, glocal = TRUE, alignment, cpu) {
    assertthat::assert_that(assertthat::is.string(cm),
                            file.exists(cm),
                            assertthat::is.flag(glocal))
    tablefile <- tempfile("cmsearch", fileext = ".dat")
    on.exit(unlink(tablefile))
    args <- c("--tblout", tablefile, "--toponly", "--noali")
    if (isTRUE(glocal)) args <- c(args, "-g")
    if (!missing(cpu)) {
        assertthat::assert_that(assertthat::is.count(cpu))
        args <- c(args, "--cpu", cpu)
    }
    if (!missing(alignment)) {
        assertthat::assert_that(assertthat::is.string(alignment))
        d <- dirname(alignment)
        if (nchar(d) > 0 && !dir.exists(d)) dir.create(d, recursive = TRUE)
        args <- c(args, "-A", alignment)
    }
    args <- c(args, cm)
    seqfile <- NULL
    if (assertthat::is.string(seq) && file.exists(seq)) {
        seqfile <- seq
    } else {
        seqfile <- tempfile("seq", fileext = ".fasta")
        if (is.character(seq)) {
            seq <- Biostrings::BStringSet(seq)
            abc <- Biostrings::uniqueLetters(seq)
            if (all(abc %in% Biostrings::DNA_ALPHABET)) {
                seq <- Biostrings::DNAStringSet(seq)
                seq <- Biostrings::RNAStringSet(seq)
            } else if (all(abc %in% Biostrings::RNA_ALPHABET)) {
                seq <- Biostrings::RNAStringSet(seq)
            } else stop("Sequence alphabet should be DNA or RNA for CMalign.")
        }
        if (methods::is(seq, "XStringSet")) {
            Biostrings::writeXStringSet(seq, seqfile)
            on.exit(unlink(seqfile))
        } else if (methods::is(seq, "ShortRead")) {
            ShortRead::writeFasta(seq, seqfile)
            on.exit(unlink(seqfile))
        } else {
            stop("'seq' should be a filename, XStringSet, ShortRead,",
                 " or character vector.")
        }
    }
    args <- c(args, seqfile)
    system2("cmsearch", args)
    readr::read_table2(
        tablefile,
        col_names = c(
            "target_name", "target_accession", "query_name", "query_accession",
            "mdl", "mdl_from", "mdl_to", "seq_from", "seq_to", "strand",
            "trunc", "pass", "gc", "bias", "score", "E_value", "inc",
            "description"),
        col_types = "ccccciiiicciddddcc", comment = "#")
}

#' Align sequences to a covariance model
#'
#' This function calls \code{cmalign} from
#' \href{http://eddylab.org/infernal/}{Infernal}.  Infernal must be installed
#' and on the path.  Not all options are included.
#'
#' One of the easiest places to obtain CMs is
#' \href{https://rfam.xfam.org/}{Rfam}.
#'
#' @param cmfile (\code{character} scalar) path to a covariance model file
#' @param seq (\code{character} scalar, names \code{character} vector,
#'        \code{\link[Biostrings]{XStringSet}}, or
#'        \code{\link[ShortRead]{ShortRead}}) sequences to align to the
#'        covariance model. This may be given as a character path to a fasta
#'        file, or the sequences as a character vector object of class
#'        \code{\link[Biostrings]{DNAStringSet}},
#'        \code{\link[Biostrings]{RNAStringSet}}, or
#'        \code{\link[ShortRead]{ShortRead}}.  For \code{cmalign}, the
#'        sequences should be known \emph{a priori} to represent the same region
#'        as the CM; to find the region in longer sequences and align it, use
#'        the \code{alignment} option of \code{\link{cmsearch}}.
#' @param glocal (\code{logical} scalar) If \code{TRUE}, align in glocal mode.
#'        See \href{http://eddylab.org/infernal/}{Infernal} documentation for
#'        more information.
#' @param cpu (\code{integer} scalar) The number of cpus to use.
#'
#' @return the aligned sequences, as returned by
#'     \code{\link{read_stockholm_msa}}.
#' @export
#'
#' @examples
#'     # align a set of unaligned 5.8S rRNA sequences to the Rfam CM.
#'     cm5_8S <- system.file(
#'         file.path("extdata", "RF00002.cm"),
#'         package = "inferrnal"
#'     )
#'     unaln <- system.file(
#'         file.path("extdata", "samp_5_8S.fasta"),
#'         package = "inferrnal"
#'     )
#'     cmalign(cm5_8S, unaln, cpu = 1)
#'     # also works if the fasta file has already been loaded
#'     unaln <- Biostrings::readRNAStringSet(unaln)
#'     cmalign(cm5_8S, unaln, cpu = 1)
cmalign <- function(cmfile, seq, glocal = TRUE, cpu) {
    assertthat::assert_that(assertthat::is.readable(cmfile),
                            assertthat::is.flag(glocal))
    args <- "cmalign"
    if (isTRUE(glocal)) args <- c(args, "-g")
    if (!missing(cpu)) {
        assertthat::assert_that(assertthat::is.count(cpu))
        args <- c(args, "--cpu", cpu)
    }

    seqfile <- NULL
    if (assertthat::is.string(seq) && file.exists(seq)) {
        seqfile <- seq
    } else {
        seqfile <- tempfile("seq", fileext = ".fasta")
        if (is.character(seq)) {
            seq <- Biostrings::BStringSet(seq)
            abc <- Biostrings::uniqueLetters(seq)
            if (all(abc %in% Biostrings::DNA_ALPHABET)) {
                seq <- Biostrings::DNAStringSet(seq)
                seq <- Biostrings::RNAStringSet(seq)
            } else if (all(abc %in% Biostrings::RNA_ALPHABET)) {
                seq <- Biostrings::RNAStringSet(seq)
            } else stop("Sequence alphabet should be DNA or RNA for CMalign.")
        }
        if (methods::is(seq, "XStringSet")) {
            Biostrings::writeXStringSet(seq, seqfile)
            on.exit(unlink(seqfile))
        } else if (methods::is(seq, "ShortRead")) {
            ShortRead::writeFasta(seq, seqfile)
            on.exit(unlink(seqfile))
        } else {
            stop("'seq' should be a filename, XStringSet, ShortRead, ",
                 "or character vector.")
        }
    }
    args <- c(args, cmfile, seqfile)
    args <- paste(args, collapse = " ")
    alnpipe <- pipe(args)
    read_stockholm_msa(alnpipe)
}
