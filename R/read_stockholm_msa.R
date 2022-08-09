# find GF (file) annotation lines and add them to the accumulator
add_gf <- function(x, acc) {
    gf <- regmatches(x, regexec("#=GF +([^ ]+) +(.+)", x))
    gf <- do.call(rbind, lapply(gf, `[`, i = 2:3))
    gf <- gf[stats::complete.cases(gf), , drop = FALSE]

    for (i in seq_len(nrow(gf))) {
        tag <- gf[i,1]
        value <- gf[i,2]
        if (is.null(acc$GF)) acc$GF <- character()
        if (tag %in% names(acc$GF)) {
            acc$GF[tag] <- paste(acc$GF[tag], value)
        } else {
            acc$GF[tag] <- value
        }
    }
    acc
}

# find GC (column) annotation lines and add them to the accumulator
add_gc <- function(x, acc) {
    gc <- regmatches(x, regexec("#=GC +([^ ]+) +(.+)", x))
    gc <- do.call(rbind, lapply(gc, `[`, i = 2:3))
    gc <- gc[stats::complete.cases(gc), , drop = FALSE]
    for (i in seq_len(nrow(gc))) {
        tag <- gc[i,1]
        value <- gc[i,2]
        if (is.null(acc$GC)) acc$GC <- character()
        if (tag %in% names(acc$GC)) {
            acc$GC[[tag]] <- paste0(acc$GC[[tag]], value)
        } else {
            acc$GC[[tag]] <- value
        }
    }
    acc
}

# find GS (sequence annotation) lines and add them to the accumulator
add_gs <- function(x, acc) {
    gs <- regmatches(x, regexec("#=GS +([^ ]+) +([^ ]+) +(.+)", x))
    gs <- do.call(rbind, lapply(gs, `[`, i = 2:4))
    gs <- gs[stats::complete.cases(gs), , drop = FALSE]
    for (i in seq_len(nrow(gs))) {
        tag <- gs[i,2]
        seq <- gs[i,1]
        value <- gs[i,3]
        if (is.null(acc$GS)) acc$GS <- list()
        if (tag %in% names(acc$GS)) {
            if (seq %in% names(acc$GS[[tag]])) {
              acc$GS[[tag]][seq] <- paste(acc$GS[[tag]][seq], value)
            } else {
              acc$GS[[tag]][seq] <- value
            }
        } else {
            acc$GS[[tag]] <- character()
            acc$GS[[tag]][seq] <- value
        }
    }
    acc
}

# find GR (residue annotation) lines and add them to the accumulator
add_gr <- function(x, acc) {

    gr <- regmatches(x, regexec("#=GR +([^ ]+) +([^ ]+) +(.+)", x))
    gr <- do.call(rbind, lapply(gr, `[`, i = 2:4))
    gr <- gr[stats::complete.cases(gr), , drop = FALSE]
    for (i in seq_len(nrow(gr))) {
        seq <- gr[i,1]
        tag <- gr[i,2]
        value <- gr[i,3]
        if (is.null(acc$GR)) acc$GR <- list()
        if (is.null(acc$GR[[tag]])) {
            acc$GR[[tag]] <- character()
            acc$GR[[tag]][seq] <- value
        } else if (seq %in% names(acc$GR[[tag]])) {
            acc$GR[[tag]][seq] <- paste0(acc$GR[[tag]][seq], value)
        } else {
            acc$GR[[tag]][seq] <- value
        }
    }
    acc
}

# find sequence lines and add them to the accumulator
add_sequences <- function(x, acc) {
    x <- regmatches(x, regexec("^(\\d+\\|)?([^#][^ ]*) +([^ ]+)$", x))
    x <- do.call(rbind, lapply(x, `[`, i = 3:4))
    x <- x[stats::complete.cases(x), , drop = FALSE]
    for (i in seq_len(nrow(x))) {
        if (x[i,1] %in% names(acc$alignment)) {
            acc$alignment[[x[i,1]]] <- paste0(acc$alignment[[x[i,1]]], x[i,2])
        } else {
            if (is.null(acc$alignment)) acc$alignment <- character()
            acc$alignment[[x[i,1]]] <- x[i,2]
        }
    }
    acc
}

parse_stockholm_msa_chunk <- function(x, pos, acc) {

    # GF (file annotation) lines
    acc <- add_gf(x, acc)

    # GC (column annotation) lines
    acc <- add_gc(x, acc)

    # GS (sequence annotation) lines
    acc <- add_gs(x, acc)

    # GR (residue annotation) lines
    acc <- add_gr(x, acc)

    # sequence lines
    acc <- add_sequences(x, acc)

    acc
}

#' Parse a Multiple Alignment in Stockholm Format
#'
#' Parses Stockholm-format multiple alignment files, including "`#=GC`"
#' lines.  Other annotations are ignored.
#'
#' @param stockholm (`character` scalar) Path to a file to parse
#' @param dna (`logical` scalar) Parse the input as DNA instead of RNA.
#'
#' @return a `list`, with elements:
#'     \describe{
#'         \item{`alignment`
#'             ([`MultipleAlignment`][Biostrings::MultipleAlignment-class])}{
#'             The alignment itself}
#'         \item{`GF` (`character`)}{Annotations which apply to the
#'             entire file.}
#'         \item{`GS` (`list` containing
#'             one [`BStringSet`][Biostrings::XStringSet-class]
#'             for each annotation)}{
#'             Annotations which apply to each sequence.}
#'         \item{`GC` (`list` of
#'             [`BString`][Biostrings::XString-class]}{
#'             Annotations which apply to each column of the alignment;
#'             in particular, "`SS_cons`" is the consensus secondary
#'             structure, and "`RF`" is the reference line.}
#'         \item{`GR` (`list` containing
#'             one [`BStringSet`][Biostrings::XStringSet-class]
#'             for each annotation)}{
#'             Annotations which apply to each residue (i.e., each column and
#'             sequence) in the alignment.}
#'     }
#' @export
#'
#' @examples
#'     msafile <- sample_rRNA_stk()
#'     msa <- read_stockholm_msa(msafile)
#'     msa$alignment
#'     # consensus secondary structure
#'     msa$GC$SS_cons
#'     # reference sequence
#'     msa$GC$RF
read_stockholm_msa <- function(stockholm, type = c("rna", "dna", "aa")) {
    assertthat::assert_that((assertthat::is.string(stockholm) &&
                                file.exists(stockholm)) ||
                                methods::is(stockholm, "connection"))
    type = match.arg(type)

    if (is.character(stockholm)) {
        stockholm <- file(stockholm)
    }

    if (!isOpen(stockholm)) {
        open(stockholm, open = "rt")
        on.exit(close(stockholm))
    }

    out <- list(alignment = list(), GF = character(), GS = list(), GR = list(),
                GC = list())
    while (TRUE) {
        lines <- readLines(stockholm, 1000, ok = TRUE)
        if (length(lines) == 0) break
        out <- parse_stockholm_msa_chunk(
            lines,
            0,
            out
        )
    }

    out$GS <- lapply(out$GS, Biostrings::BStringSet)
    out$GR <- lapply(out$GR, Biostrings::BStringSet)
    out$GC <- lapply(out$GC, Biostrings::BString)

    if (type == "dna") {
        out$alignment <- Biostrings::DNAMultipleAlignment(unlist(out$alignment))
    } else if (type == "rna") {
        out$alignment <- Biostrings::RNAMultipleAlignment(unlist(out$alignment))
    } else if (type == "aa") {
        out$alignment <- Biostrings::AAMultipleAlignment(unlist(out$alignment))
    } else {
      stop("unknown sequence type: ", type)
    }
    out
}
