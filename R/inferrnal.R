#' Search for Covariance Models (CM) in a set of sequences.
#'
#' This function calls "\code{cmsearch}" from
#' \href{http://eddylab.org/infernal/}{Infernal}.  Infernal must be installed.
#' Many parameters are not included (yet!), and the function is focused on
#' retrieving the hits table and, optionally, producing an alignment.
#'
#' @param cm (character of length 1) Path to the covariance model (.cm) file.
#'     The covariance model must include calibration data from running
#'     "\code{cmcalibrate}".
#' @param seq (filename, character vector,
#'     \code{\link[Biostrings]{XStringSet-class}}, or
#'     \code{\link[ShortRead]{ShortRead-class}}) Sequences to search with the CM.
#'     If a filename, the file can be of any type supported by Infernal.
#' @param glocal (logical of length 1) Whether to run the search in glocal mode
#'     (global with respect to the model, local with respect to the sequence).
#'     When \code{TRUE}, the search is faster, but will fail to find matches
#'     with only partially homologous structure.
#' @param Z (\code{numeric} scalar) Effective search space size in Mb for the
#'     purposes of calculating E-values.
#' @param output (\code{character} filename) File to send the human-readable
#'     output to.
#' @param alignment (filename) A file to save the aligned hits to.  If given,
#'     the alignment is saved in Stockholm format with annotations for secondary
#'     structure, posterior probablility, etc.
#' @param acc (\code{logical} scalar) Use accessions instead of names in the
#'     main output.
#' @param noali (\code{logical} scalar) Omit the alignment section from the main
#'     output.
#' @param notextw (\code{logical} scalar) Unlimit the length of each line in the
#'    main output.
#' @param textw (\code{numeric} scalar) Set the main output’s line length limit
#'     in characters per line. The default is 120.
#' @param verbose (\code{logical} scalar)
#' @param E (\code{numeric} scalar) Maximum E-value for reporting in per-target
#'     output. Default: 10.0
#' @param T (\code{numeric} scalar) Maximum bit score for reporting in
#'     per-target reporting. Default: NULL
#' @param incE (\code{numeric} scalar) Maximum E-value for hit inclusion.
#'     Default: 0.01
#' @param incT (\code{numeric} scalar) Maximum bit score for hit inclusion.
#'     Default: NULL
#' @param cut_ga (\code{logical} scalar) Use the GA (gathering) bit scores
#'     defined in the CM to set hit reporting and inclusion thresholds.
#' @param cut_nc (\code{logical} scalar) Use the NC (noise cutoff) bit score
#'     thresholds defined in the CM to set hit reporting and inclusion
#'     thresholds.
#' @param cut_tc (\code{logical} scalar) Use the TC (trusted cutoff) bit score
#'   thresholds defined in the CM to set hit reporting and inclusion thresholds.
#' @param filter_strategy (\code{character} string) Filtering strategy for the
#'     acceleration pipeline. Options, from slowest (most sensitive) to fastest
#'     (least sensitive) are "max", "nohmm", "mid", "default", "rfam",
#'     "hmmonly".
#' @param FZ (\code{numeric} scalar) Effective database size in Mb for
#'     determining filter thresholds.
#' @param Fmid (\code{numeric} scalar) HMM filter thresholds for the "mid"
#'     filtering strategy.  Default 0.02
#' @param notrunc (\code{logical} scalar) Turn off truncated hit detection.
#' @param anytrunc (\code{logical} scalar) Allow truncated hits to begin and end
#'     at any position in a target sequence.
#' @param nonull3 (\code{logical} scalar) Turn off the null3 CM score
#'     corrections for biased composition.
#' @param mxsize (\code{numeric} scalar) Maximum allowable CM DP matrix size in
#'   megabytes. Default: 128
#' @param smxsize (\code{numeric} scalar) Maximum allowable CM search DP matrix
#'     size in megabytes. Default: 128
#' @param cyk (\code{logical} scalar) Use the CYK algorithm, not Inside, to
#'     determine the final score of all hits.
#' @param acyk (\code{logical} scalar) Use the CYK algorithm to align hits.
#' @param wcx (\code{numeric} scalar) Expected maximum length of a hit, as a
#'     multiplier of consensus length of the model.
#' @param toponly (\code{logical} scalar) Only search the top (Watson) strand of
#'     target sequences.
#' @param bottomonly (\code{logical} scalar) Only search the bottom (Crick)
#'     strand of target sequences.
#' @param tformat (\code{logical} scalar) Format of the sequences in "seq".
#'     Options are "fasta", "embl", "genbank", "ddbj", "stockholm", "pfam",
#'     "a2m", "afa", "clustal", and "phylip". Default: autodetect
#' @param cpu (\code{integer} scalar) The number of threads to use in the
#'     search.
#' @param stall (\code{logical} scalar) For debugging the MPI master/worker
#'     version: pause after start, to enable the developer toattach debuggers to
#'     the running master and worker(s) processes
#' @param mpi (\code{logical} scalar) Run  in  MPI  master/worker  mode, using
#'     mpirun.
#' @param quiet (\code{logical} scalar) Suppress standard output of `cmsearch`,
#'     which can be long.
#'
#' @return a \code{\link{data.frame}} with columns:
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
cmsearch <- function(
    cm,
    seq,
    glocal = TRUE,
    Z = NULL,
    output = NULL,
    alignment = NULL,
    acc = FALSE,
    noali = TRUE,
    notextw = FALSE,
    textw = NULL,
    verbose = FALSE,
    E = NULL,
    T = NULL,
    incE = NULL,
    incT = NULL,
    cut_ga = FALSE,
    cut_nc = FALSE,
    cut_tc = FALSE,
    filter_strategy = NULL,
    FZ = NULL,
    Fmid = NULL,
    notrunc = FALSE,
    anytrunc = FALSE,
    nonull3 = FALSE,
    mxsize = NULL,
    smxsize = NULL,
    cyk = FALSE,
    acyk = FALSE,
    wcx = NULL,
    toponly = TRUE,
    bottomonly = FALSE,
    tformat = NULL,
    cpu = NULL,
    stall = FALSE,
    mpi = FALSE,
    quiet = TRUE) {
    assertthat::assert_that(assertthat::is.string(cm),
                            file.exists(cm),
                            assertthat::is.flag(glocal))
    tablefile <- tempfile("cmsearch", fileext = ".dat")
    on.exit(unlink(tablefile))
    args <- c(
        flag_opt(glocal, "g"),
        nonneg_float_opt(Z),
        string_opt(output, "o"),
        string_opt(alignment, "A"),
        string_opt(tablefile, "tblout"),
        flag_opt(acc),
        flag_opt(noali),
        flag_opt(notextw),
        count_opt(textw),
        flag_opt(verbose),
        nonneg_float_opt(E),
        nonneg_float_opt(T),
        nonneg_float_opt(incE),
        nonneg_float_opt(incT),
        flag_opt(cut_ga),
        flag_opt(cut_nc),
        flag_opt(cut_tc),
        multiflag_opt(
            filter_strategy,
            choices = c("max", "nohmm", "mid", "default", "rfam", "hmmonly")
        ),
        nonneg_float_opt(FZ),
        fraction_opt(Fmid),
        flag_opt(notrunc),
        flag_opt(anytrunc),
        flag_opt(nonull3),
        nonneg_float_opt(mxsize),
        nonneg_float_opt(smxsize),
        flag_opt(cyk),
        flag_opt(acyk),
        nonneg_float_opt(wcx),
        flag_opt(toponly),
        flag_opt(bottomonly),
        string_opt(
            tformat,
            choices = c(
                "fasta", "embl", "genbank", "ddbj", "stockholm", "pfam", "a2m",
                "afa", "clustal", "phylip"
            )
        ),
        count_opt(cpu),
        flag_opt(stall),
        flag_opt(mpi),
        cm
    )

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
            stop(
                "'seq' should be a filename, XStringSet, ShortRead,",
                " or character vector."
            )
        }
    }
    args <- c(args, seqfile)
    system2(
        "cmsearch",
        args,
        stdout = if (isTRUE(quiet)) FALSE else "",
        stderr = if (isTRUE(quiet)) FALSE else ""
    )
    utils::read.table(
        tablefile,
        col.names = c(

            #character
            "target_name",
            "target_accession",
            "query_name",
            "query_accession",
            "mdl",

            #integer
            "mdl_from",
            "mdl_to",
            "seq_from",
            "seq_to",

            # character
            "strand",
            "trunc",

            # integer
            "pass",

            # numeric
            "gc",
            "bias",
            "score",
            "E_value",

            #character
            "inc",
            "description"
        ),
        colClasses = c(
            rep("character", 5),
            rep("integer", 4),
            rep("character", 2),
            "integer",
            rep("numeric", 4),
            rep("character", 2)
        ),
        comment.char = "#")
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
#' @param cmfile (\code{character} filename) path to a covariance model file
#' @param seq (\code{character} filename, \code{character} vector,
#'     \code{\link[Biostrings]{XStringSet-class}}, or
#'     \code{\link[ShortRead]{ShortRead-class}}) sequences to align to the
#'     covariance model. This may be given as a character path to a fasta
#'     file, or the sequences as a character vector object of class
#'     \code{\link[Biostrings]{DNAStringSet-class}},
#'     \code{\link[Biostrings]{RNAStringSet-class}}, or
#'     \code{\link[ShortRead]{ShortRead-class}}.  For \code{cmalign}, the
#'     sequences should be known \emph{a priori} to represent the same region
#'     as the CM; to find the region in longer sequences and align it, use
#'     the \code{alignment} option of \code{\link{cmsearch}}.
#' @param global (\code{logical} scalar) If \code{TRUE}, align in global mode.
#'     See \href{http://eddylab.org/infernal/}{Infernal} documentation for
#'     more information.
#' @param algorithm (\code{character} string) Alignment algorithm. Options are
#'     "optacc" and "cyk". Default: "optacc"
#' @param sample (\code{logical} scalar) Sample an alignment from the posterior
#'     distribution of alignments.
#' @param seed (\code{integer} scalar) Random seed for sampling an alignment.
#' @param notrunc (\code{logical} scalar) Turn off truncated alignment
#'     algorithms.
#' @param sub (\code{logical} scalar) urn on the sub model construction and
#'     alignment procedure
#' @param hbanded (\code{logical} scalar) Accelerate alignment by pruning away
#'     regions of the CMDP matrix that are deemed negligible by an HMM.
#'     Default: TRUE
#' @param tau (\code{numeric} scalar) Tail loss probability used during HMM band
#'     calculation. Default: 1E-7
#' @param mxsize (\code{numeric} scalar) The maximum DP matrix size, in Mb.
#'     Maximum potential memory usage is approximately cpu*mxsize, although
#'     this is usually not realized.
#'     See \href{http://eddylab.org/infernal/}{Infernal} documentation for
#'     more information.
#' @param fixedtau (\code{logical} scalar) Turn off HMM band tightening.
#' @param maxtau (\code{numeric} scalar) Maximum allowed value for tau
#'     during band tightening. Default: 0.05
#' @param small (\code{logical} scalar) Use the divide and conquer CYK alignment
#'     algorithm, greatly reducing memory consumption.
#' @param sfile (\code{character} filename) Dump per-sequence alignment score
#'     and timing information to file.
#' @param tfile (\code{character} filename) Dump tabular sequence tracebacks for
#'     each individual sequence to a file.
#' @param ifile (\code{character} filename) Dump per-sequence insert information
#'     to file.
#' @param elfile (\code{character} filename) Dump per-sequence EL state (local
#'     end) insert information to file.
#' @param mapali (\code{character} filename) Read the alignment from the file
#'     used
#'     to build the model aligns it as a single object to the CM, along with
#'     sequences in \code{"seq"}.
#' @param mapstr (\code{logical} scalar) Must be used in combination with
#'     \code{mapali}. Propogate structural information for any pseudoknots
#'     that exist in \code{mapali} to the output alignment.
#' @param dnaout (\code{logical} scalar) Output the alignments as DNA sequence
#'     alignments, instead of RNA ones.
#' @param noprob (\code{logical} scalar) Do not annotate the output alignment
#'     with posterior probabilities.
#' @param matchonly (\code{logical} scalar) Only include match columns in the
#'     output alignment, do not include any insertions relativeto the
#'     consensus model.
#' @param ileaved (\code{logical} scalar) Output the alignment in interleaved
#'     Stockholm format of a fixed width that may be more convenient for
#'     examination.
#' @param regress (\code{character} filename) Save an additional copy of the
#'     output alignment with no author information to file.
#' @param verbose (\code{logical} scalar) Output additional information in the
#'     tabular scores output.
#' @param cpu (\code{integer} scalar) The number of cpus to use.
#' @param mpi (\code{logical} scalar) Run as an MPI parallel program.
#' @param extra (\code{character} vector) Additional advanced options to pass
#'     to `cmalign`.
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
cmalign <- function(
    cmfile,
    seq,
    global = TRUE,
    algorithm = NULL,
    sample = FALSE,
    seed = NULL,
    notrunc = FALSE,
    sub = FALSE,
    hbanded = TRUE,
    tau = NULL,
    mxsize = NULL,
    fixedtau = FALSE,
    maxtau = NULL,
    small = FALSE,
    sfile = NULL,
    tfile = NULL,
    ifile = NULL,
    elfile = NULL,
    mapali = NULL,
    mapstr = FALSE,
    dnaout = FALSE,
    noprob = FALSE,
    matchonly = FALSE,
    ileaved = FALSE,
    regress = NULL,
    verbose = FALSE,
    cpu = NULL,
    mpi = FALSE,
    extra = NULL,
    glocal = global
) {
    assertthat::assert_that(assertthat::is.readable(cmfile))
    if (!missing(glocal))  {
        warning("Use of 'glocal' in cmalign is deprecated. ",
                "Please use 'global' instead.")
        if (!missing(global)) {
            stop(
                "'global' and 'glocal' should not both be specified in",
                " cmalign."
            )
        }
        global <- glocal
    }

    args <- c(
        "cmalign",
        flag_opt(global, "g"),
        multiflag_opt(algorithm, c("optacc", "cyk")),
        flag_opt(sample),
        count_opt(seed),
        flag_opt(notrunc),
        flag_opt(sub),
        flag_opt(hbanded, "nonbanded", invert = TRUE),
        fraction_opt(tau),
        count_opt(mxsize),
        flag_opt(fixedtau),
        fraction_opt(maxtau),
        flag_opt(small),
        string_opt(sfile),
        string_opt(tfile),
        string_opt(ifile),
        string_opt(elfile),
        infile_opt(mapali, ),
        flag_opt(mapstr),
        flag_opt(dnaout),
        flag_opt(noprob),
        flag_opt(matchonly),
        flag_opt(ileaved),
        string_opt(regress),
        flag_opt(verbose),
        count_opt(cpu),
        flag_opt(mpi),
        extra
    )

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
            stop(
                "'seq' should be a filename, XStringSet, ShortRead, ",
                "or character vector."
            )
        }
    }
    args <- c(args, cmfile, seqfile)
    args <- paste(shQuote(args), collapse = " ")
    alnpipe <- pipe(args)
    read_stockholm_msa(alnpipe, dna = dnaout)
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
#' @param name (\code{character} string) name of the CM (option \code{-N} to
#'     cmbuild). The default uses the
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
#' @param weights (\code{character}; one of \code{"wpb"}, \code{"wgsc"},
#'     \code{"wnone"}, \code{"wgiven"}, or \code{"wblosum"}) sequence weighting
#'     method (options \code{--wpb}, \code{--wgsc},
#'     \code{--wnone}, \code{--wgiven}, and \code{--wblosum} to cmbuild).
#' @param wid (\code{numeric} scalar) percent identity for clustering when
#'     \code{weights = "blosum"}.
#' @param eff_num (\code{character}; one of \code{"eent"} or \code{"enone"})
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
#'     sampling.
#' @param cyk (\code{logical} scalar) use the CYK alignment algorithm with
#'     \code{refine}.
#' @param notrunc (\code{logical} scalar) turn off the truncated alignment
#'     algorithm with \code{refine}.
#' @param extra (\code{character}) additional arguments to pass to cmbuild.
#' @param quiet (\code{logical} scalar) do not print cmbuild output to the
#'     console.
#'
#' @return \code{NULL}, invisibly
#' @export
#'
#' @examples
#'     cmbuild(
#'         msafile = stk_5_8S(),
#'         cmfile_out = "/dev/null",
#'         force = TRUE,
#'         quiet = FALSE
#'     )
#'
cmbuild <- function(
    msafile,
    cmfile_out,
    name = NULL,
    force = FALSE,
    summary_file = NULL,
    reannotated_msa = NULL,
    consensus_method = NULL,
    symfrac = NULL,
    rsearch = NULL,
    null = NULL,
    prior = NULL,
    weights = NULL,
    wid = NULL,
    eff_num = NULL,
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
    quiet = TRUE
) {
    args <- c(
        string_opt(name, "n"),
        flag_opt(force, "F"),
        string_opt(summary_file, "O"),
        multiflag_opt(consensus_method, choices = c("fast", "hand", "noss")),
        fraction_opt(symfrac),
        infile_opt(rsearch),
        infile_opt(null),
        infile_opt(prior),
        multiflag_opt(weights, choices = c("wpb", "wgsc", "wgiven", "wblosum")),
        percent_opt(wid),
        multiflag_opt(eff_num, choices = c("eent", "enone")),
        nonneg_float_opt(ere),
        nonneg_float_opt(eminseq),
        nonneg_float_opt(emaxseq),
        nonneg_float_opt(ehmmre),
        nonneg_float_opt(eset),
        nonneg_float_opt(p7ere),
        flag_opt(p7ml),
        string_opt(refine),
        flag_opt(local),
        flag_opt(gibbs),
        count_opt(seed),
        flag_opt(cyk),
        flag_opt(notrunc),
        extra,
        cmfile_out,
        msafile
    )

    assertthat::assert_that(assertthat::is.flag(quiet))

    return <- system2(
        command = "cmbuild",
        args = args,
        stdout = if (isTRUE(quiet)) FALSE else "",
        stderr = if (isTRUE(quiet)) FALSE else ""
    )

    if (return != 0) stop("cmbuild failed with exit code ", return)
    invisible(NULL)
}
