
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
#'     msafile <- sample_rRNA_stk()
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
#'     cm <- cm_5_8S()
#'     sampfasta <- sample_rRNA_fasta()
#'     cmsearch(cm = cm, seq = sampfasta, cpu = 1)
#'     # also works if the fasta file has already been loaded
#'     samp <- Biostrings::readDNAStringSet(sampfasta)
#'     cmsearch(cm = cm, seq = samp, cpu = 1)
cmsearch <- function(cm, seq, glocal = TRUE, alignment = NULL, cpu = NULL) {
    assertthat::assert_that(assertthat::is.string(cm),
                            file.exists(cm),
                            assertthat::is.flag(glocal))
    tablefile <- tempfile("cmsearch", fileext = ".dat")
    on.exit(unlink(tablefile))
    args <- c("--tblout", tablefile, "--toponly", "--noali")
    if (isTRUE(glocal)) args <- c(args, "-g")
    if (!is.null(cpu)) {
        assertthat::assert_that(assertthat::is.count(cpu))
        args <- c(args, "--cpu", cpu)
    }
    if (!is.null(alignment)) {
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
#' @param global (\code{logical} scalar) If \code{TRUE}, align in global mode.
#'        See \href{http://eddylab.org/infernal/}{Infernal} documentation for
#'        more information.
#' @param cpu (\code{integer} scalar) The number of cpus to use.
#' @param mxsize (\code{double} scalar) The maximum DP matrix size, in Mb.
#'        Maximum potential memory usage is approximately cpu*mxsize, although
#'        this is usually not realized.
#'        See \href{http://eddylab.org/infernal/}{Infernal} documentation for
#'        more information.
#' @param glocal (\code{logical} scalar) (Deprecated) see "global".
#'
#' @return the aligned sequences, as returned by
#'     \code{\link{read_stockholm_msa}}.
#' @export
#'
#' @examples
#'     # align a set of unaligned 5.8S rRNA sequences to the Rfam CM.
#'     cm <- cm_5_8S()
#'     unaln <- sample_rRNA_5_8S()
#'     cmalign(cm, unaln, cpu = 1)
#'     # also works if the fasta file has already been loaded
#'     unaln <- Biostrings::readRNAStringSet(unaln)
#'     cmalign(cm, unaln, cpu = 1)
cmalign <- function(cmfile, seq, global = TRUE, cpu = NULL, mxsize = NULL,
                    glocal = global) {
    assertthat::assert_that(assertthat::is.readable(cmfile),
                            assertthat::is.flag(glocal))
    args <- "cmalign"
    if (!missing(glocal))  {
        warning("Use of 'glocal' in cmalign is deprecated. ",
                "Please use 'global' instead.")
        if (!missing(global)) {
            stop("'global' and 'glocal' should not both be specified in",
                 " cmalign.")
        }
        global <- glocal
    }
    if (isTRUE(global)) args <- c(args, "-g")
    if (!is.null(cpu)) {
        assertthat::assert_that(assertthat::is.count(cpu))
        args <- c(args, "--cpu", cpu)
    }
    if (!is.null(mxsize)) {
        assertthat::assert_that(
            assertthat::is.number(mxsize),
            mxsize > 0)
        args <- c(args, "--mxsize", mxsize)
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

#' Construct covariance model(s) from structually annotated alignment(s)
#'
#' Calls the standalone program \code{cmbuild} from the
#' \href{http://eddylab.org/infernal/}{Infernal} package, which must be
#' installed.  For more information about options, see the
#' \href{http://eddylab.org/infernal/Userguide.pdf}{Infernal documentation}.
#'
#' @param msafile (\code{character} filename) filename of MSA file to read.
#' @param cmfile_out (\code{character} filename) filename to write CM to.
#' @param name (\code{character} string) name of the CM (option \code{-N} to cmbuild). The default uses the
#'     name(s) given in the alignment file, or the name of the alignment file.
#' @param force (\code{logical} scalar) overwrite \code{cmfile_out} if it
#'     exists
#' @param summary_file (\code{character} filename) filename to write summary
#'     output to.
#' @param reannotated_msa (\code{character} filename) filename to write a
#'     reannotated alignment(s) to (option \code{-O} to cmbuild)
#' @param consensus_method (one of \code{"fast"}, \code{"hand"}, or
#'     \code{"noss"}) method to define consensus columns in alignment. (options
#'     \code{--fast}, \code{--hand}, and \code{--noss} to cmbuild)
#' @param symfrac (\code{integer} scalar) residue fraction threshold necessary
#'     to define a consensus column.
#' @param rsearch (\code{character} filename) RIBOSUM matrix to use in
#'     parameterizing emission scores.
#' @param null (\code{character} filename) null model file.
#' @param prior (\code{character} filename) Dirichlet prior file.
#' @param weights (\code{character}; one of \code{"pb"}, \code{"gsc"},
#'     \code{"none"}, \code{"given"}, or \code{"blosum"}) sequence weighting
#'     method (options \code{--wpb}, \code{--wgsc},
#'     \code{--wnone}, \code{--wgiven}, and \code{wblosum} to cmbuild).
#' @param wid (\code{numeric} scalar) percent identity for clustering when
#'     \code{weights = "blosum"}.
#' @param eff_num (\code{character}; one of \code{"ent"} or \code{"none"})
#'     entropy weighting strategy (options \code{--eent} or \code{--enone} to
#'     cmbuild).
#' @param ere (\code{numeric} scalar) target mean match state relative entropy.
#' @param eminseq (\code{numeric} scalar) minimum allowed effective sequence
#'     number.
#' @param emaxseq (\code{numeric} scalar) maximum allowed effective sequence
#'     number.
#' @param ehmmre (\code{numeric} scalar) target HMM mean match state relative
#'     entropy.
#' @param eset (\code{numeric} scalar) effective sequence number for entropy
#'     weighting.
#' @param p7ere (\code{numeric} scalar) target mean match state relative entropy
#'     for the filter p7 HMM.
#' @param p7ml (\code{logical} scalar) use a maximum liklihood p7 HMM built from
#'     the CM.
#' @param refine (\code{character} filename) if given, the alignment is
#'     iteratively refined, and the final alignment is written to this file.
#' @param local (\code{logical} scalar) use local alignment with \code{refine}
#'     (option \code{-l} to cmbuild).
#' @param gibbs (\code{logical} scalar) use Gibbs sampling instead of
#'     maximum likelihood with \code{refine}.
#' @param seed (\code{integer} scalar) pseudorandom number seed for Gibbs
#'    sampling.
#' @param cyk (\code{logical} scalar) use the CYK alignment algorithm with
#'    \code{refine}.
#' @param notrunc (\code{logical} scalar) turn off the truncated alignment
#'    algorithm with \code{refine}.
#' @param extra (\code{character}) additional arguments to pass to cmbuild.
#' @param verbose (\code{logical} scalar) print cmbuild output to the console.
#'
#' @return \code{NULL}, invisibly
#' @export
#'
#' @examples
#'
#' # requires that LSUx is installed
#'
#'
#'
cmbuild <- function(
    msafile,
    cmfile_out,
    name = NULL,
    force = FALSE,
    summary_file = NULL,
    reannotated_msa = NULL,
    consensus_method = c("fast", "hand", "noss"),
    symfrac = NULL,
    rsearch = NULL,
    null = NULL,
    prior = NULL,
    weights = c("pb", "gsc", "given", "blosum"),
    wid = NULL,
    eff_num = c("ent", "none"),
    ere = NULL,
    eminseq = NULL,
    emaxseq = NULL,
    ehmmre = NULL,
    eset = NULL,
    p7ere = NULL,
    p7ml = FALSE,
    refine = NULL,
    local = FALSE,
    gibbs = FALSE,
    seed = NULL,
    cyk = FALSE,
    notrunc = FALSE,
    extra = NULL,
    verbose = FALSE
) {
    args <- character()

    if (!is.null(name)) {
        assertthat::assert_that(assertthat::is.string(name))
        args <- c(args, "-n", name)
    }

    assertthat::assert_that(assertthat::is.flag(force))
    if (isTRUE(force)) {
        args <- c(args, "-F")
    }

    if (!is.null(summary_file)) {
        assertthat::assert_that(assertthat::is.string(summary_file))
        args <- c(args, "-O", summary_file)
    }

    if (!missing(consensus_method)) {
        consensus_method = match.arg(consensus_method)
        args <- c(args, paste0("--", consensus_method))
    }

    if (!is.null(symfrac)) {
        assertthat::assert_that(
            assertthat::is.number(symfrac),
            symfrac >= 0,
            symfrac <= 1
        )
        args <- c(args, "--symfrac", symfrac)
    }

    if (!is.null(rsearch)) {
        assertthat::assert_that(
            assertthat::is.string(rsearch),
            assertthat::is.readable(rsearch)
        )
        args <- c(args, "--rsearch", rsearch)
    }

    if (!is.null(null)) {
        assertthat::assert_that(
            assertthat::is.string(null),
            assertthat::is.readable(null)
        )
        args <- c(args, "--null", null)
    }

    if (!is.null(prior)) {
        assertthat::assert_that(
            assertthat::is.string(prior),
            assertthat::is.readable(prior)
        )
        args <- c(args, "--prior", prior)
    }

    if (!missing(weights)) {
        weights = match.arg(weights)
        args <- c(args, paste0("--w", weights))
    }

    if (!is.null(wid)) {
        assertthat::assert_that(
            assertthat::is.number(wid),
            wid >= 0,
            wid <= 100
        )
        args <- c(args, "--wid", wid)
    }

    if (!missing(eff_num)) {
        eff_num = match.arg(eff_num)
        args <- c(args, paste0("--e", eff_num))
    }

    if (!is.null(ere)) {
        assertthat::assert_that(
            assertthat::is.number(ere),
            ere >= 0
        )
        args <- c(args, "--ere", ere)
    }

    if (!is.null(eminseq)) {
        assertthat::assert_that(
            assertthat::is.number(eminseq),
            eminseq >= 0
        )
        args <- c(args, "--eminseq", eminseq)
    }

    if (!is.null(emaxseq)) {
        assertthat::assert_that(
            assertthat::is.number(emaxseq),
            emaxseq >= 0
        )
        args <- c(args, "--emaxseq", emaxseq)
    }

    if (!is.null(ehmmre)) {
        assertthat::assert_that(
            assertthat::is.number(ehmmre),
            ehmmre >= 0
        )
        args <- c(args, "--ehmmre", ehmmre)
    }

    if (!is.null(eset)) {
        assertthat::assert_that(
            assertthat::is.number(eset),
            eset >= 0
        )
        args <- c(args, "--eset", eset)
    }

    if (!is.null(p7ere)) {
        assertthat::assert_that(
            assertthat::is.number(p7ere),
            p7ere >= 0
        )
        args <- c(args, "--p7ere", p7ere)
    }

    assertthat::assert_that(assertthat::is.flag(p7ml))
    if (isTRUE(p7ml)) {
        args <- c(args, "--p7ml")
    }

    if (!is.null(refine)) {
        assertthat::assert_that(assertthat::is.string(refine))
        args <- c(args, "--refine", refine)
    }

    assertthat::assert_that(assertthat::is.flag(local))
    if (isTRUE(local)) {
        args <- c(args, "-l")
    }

    assertthat::assert_that(assertthat::is.flag(gibbs))
    if (isTRUE(gibbs)) {
        args <- c(args, "--gibbs")
    }

    if (!is.null(seed)) {
        assertthat::assert_that(
            assertthat::is.count(seed)
        )
        args <- c(args, "--seed", seed)
    }

    assertthat::assert_that(assertthat::is.flag(cyk))
    if (isTRUE(cyk)) {
        args <- c(args, "--cyk")
    }

    assertthat::assert_that(assertthat::is.flag(notrunc))
    if (isTRUE(notrunc)) {
        args <- c(args, "--notrunc")
    }

    if (!is.null(extra)) {
        assertthat::assert_that(is.character(extra))
        args <- c(args, extra)
    }

    args <- c(args, cmfile_out, msafile)

    assertthat::assert_that(assertthat::is.flag(verbose))

    return <- system2(
        command = "cmbuild",
        args = args,
        stdout = if (verbose) "" else FALSE,
        stderr = if (verbose) "" else FALSE
    )

    if (return != 0) stop("cmbuild failed with exit code ", return)
    invisible(NULL)
}
