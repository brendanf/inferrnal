#' Search for Covariance Models (CM) in a set of sequences.
#'
#' This function calls "`cmsearch`" from
#' [Infernal](http://eddylab.org/infernal/).  Infernal must be installed.
#' The function is focused on
#' retrieving the hits table and, optionally, producing an alignment.
#'
#' @param cm (character of length 1) Path to the covariance model (.cm) file.
#'     The covariance model must include calibration data from running
#'     "`cmcalibrate`".
#' @param seq (filename, character vector,
#'     [`Biostrings::XStringSet-class`], or
#'     [`ShortRead::ShortRead-class`] Sequences to search with the CM.
#'     If a filename, the file can be of any type supported by Infernal.
#' @param glocal (`logical` of length 1) Whether to run the search in glocal
#'     mode
#'     (global with respect to the model, local with respect to the sequence).
#'     When `TRUE`, the search is faster, but will fail to find matches
#'     with only partially homologous structure.
#' @param Z (`numeric` scalar) Effective search space size in Mb for the
#'     purposes of calculating E-values.
#' @param output (`character` filename) File to send the human-readable
#'     output to.
#' @param alignment (`character` filename) A file to save the aligned hits to.
#'     If given,
#'     the alignment is saved in Stockholm format with annotations for secondary
#'     structure, posterior probablility, etc.
#' @param acc (`logical` scalar) Use accessions instead of names in the
#'     main output.
#' @param noali (`logical` scalar) Omit the alignment section from the main
#'     output.
#' @param notextw (`logical` scalar) Unlimit the length of each line in the
#'    main output.
#' @param textw (`numeric` scalar) Set the main output’s line length limit
#'     in characters per line. Default: `120`.
#' @param verbose (`logical` scalar)
#' @param E (`numeric` scalar) Maximum E-value for reporting in per-target
#'     output. Default: `10.0`
#' @param t (`numeric` scalar) Maximum bit score for reporting in
#'     per-target reporting. Default: `NULL`
#' @param incE (`numeric` scalar) Maximum E-value for hit inclusion.
#'     Default: `0.01`
#' @param incT (`numeric` scalar) Maximum bit score for hit inclusion.
#'     Default: `NULL`
#' @param cut_ga (`logical` scalar) Use the GA (gathering) bit scores
#'     defined in the CM to set hit reporting and inclusion thresholds.
#' @param cut_nc (`logical` scalar) Use the NC (noise cutoff) bit score
#'     thresholds defined in the CM to set hit reporting and inclusion
#'     thresholds.
#' @param cut_tc (`logical` scalar) Use the TC (trusted cutoff) bit score
#'   thresholds defined in the CM to set hit reporting and inclusion thresholds.
#' @param filter_strategy (`character` string) Filtering strategy for the
#'     acceleration pipeline. Options, from slowest (most sensitive) to fastest
#'     (least sensitive) are `"max"`, `"nohmm"`, `"mid"`, `"default"`, `"rfam"`,
#'     `"hmmonly"`.
#' @param FZ (`numeric` scalar) Effective database size in Mb for
#'     determining filter thresholds.
#' @param Fmid (`numeric` scalar) HMM filter thresholds for the "mid"
#'     filtering strategy.  Default: `0.02`
#' @param notrunc (`logical` scalar) Turn off truncated hit detection.
#' @param anytrunc (`logical` scalar) Allow truncated hits to begin and end
#'     at any position in a target sequence.
#' @param nonull3 (`logical` scalar) Turn off the null3 CM score
#'     corrections for biased composition.
#' @param mxsize (`numeric` scalar) Maximum allowable CM DP matrix size in
#'   megabytes. Default: `128`
#' @param smxsize (`numeric` scalar) Maximum allowable CM search DP matrix
#'     size in megabytes. Default: `128`
#' @param cyk (`logical` scalar) Use the CYK algorithm, not Inside, to
#'     determine the final score of all hits.
#' @param acyk (`logical` scalar) Use the CYK algorithm to align hits.
#' @param wcx (`numeric` scalar) Expected maximum length of a hit, as a
#'     multiplier of consensus length of the model.
#' @param toponly (`logical` scalar) Only search the top (Watson) strand of
#'     target sequences.
#' @param bottomonly (`logical` scalar) Only search the bottom (Crick)
#'     strand of target sequences.
#' @param tformat (`character` string) Format of the sequences in "`seq`".
#'     Options are `"fasta"`, `"embl"`, `"genbank"`, `"ddbj"`, `"stockholm"`,
#'     `"pfam"`, `"a2m"`, `"afa"`, `"clustal"`, and `"phylip"`.
#'     Default: autodetect
#' @param cpu (`integer` scalar) The number of threads to use in the
#'     search.
#' @param stall (`logical` scalar) For debugging the MPI master/worker
#'     version: pause after start, to enable the developer toattach debuggers to
#'     the running master and worker(s) processes
#' @param mpi (`logical` scalar) Run  in  MPI  master/worker mode, using
#'     mpirun.
#' @param quiet (`logical` scalar) Suppress standard output of `cmsearch`,
#'     which can be long.
#'
#' @return a [`base::data.frame`] with columns:
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
    t = NULL,
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
        nonneg_float_opt(t, "T"),
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
            if (!requireNamespace("ShortRead", quietly = TRUE)) {
                stop("'seq' is a ShortRead object but the ShortRead package",
                    "is unavailable")
            }
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
    header <- readLines(tablefile, 2)
    gaps <- gregexpr(" ", header[2])[[1]]
    w <- c(gaps, 500) - c(0, gaps)
    utils::read.fwf(
        tablefile,
        widths = w,
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
        comment.char = "#",
        strip.white = TRUE
    )
}


#' Align sequences to a covariance model
#'
#' This function calls `cmalign` from
#'[Infernal](http://eddylab.org/infernal/).  Infernal must be installed
#' and on the path.  Not all options are included.
#'
#' One of the easiest places to obtain CMs is
#' [Rfam](https://rfam.xfam.org/).
#'
#' @param cmfile (`character` filename) path to a covariance model file
#' @param seq (`character` filename, `character` vector,
#'     [`Biostrings::XStringSet-class`], or
#'     [`ShortRead::ShortRead-class`] sequences to align to the
#'     covariance model. This may be given as a character path to a fasta
#'     file, the sequences as a character vector, or an object of class
#'     [`DNAStringSet`][Biostrings::XStringSet-class],
#'     [`RNAStringSet`][Biostrings::XStringSet-class], or
#'     [`ShortRead`][ShortRead::ShortRead-class].  For `cmalign`, the
#'     sequences should be known _a priori_ to represent the same region
#'     as the CM; to find the region in longer sequences and align it, use
#'     the `alignment` option of [cmsearch()].
#' @param global (`logical` scalar) If `TRUE`, align in global mode.
#'     See [Infernal](http://eddylab.org/infernal/) documentation for
#'     more information.
#' @param algorithm (`character` string) Alignment algorithm. Options are
#'     `"optacc"` and `"cyk"`. Default: `"optacc"`
#' @param sample (`logical` scalar) Sample an alignment from the posterior
#'     distribution of alignments.
#' @param seed (`integer` scalar) Random seed for sampling an alignment.
#' @param notrunc (`logical` scalar) Turn off truncated alignment
#'     algorithms.
#' @param sub (`logical` scalar) turn on the sub model construction and
#'     alignment procedure
#' @param hbanded (`logical` scalar) Accelerate alignment by pruning away
#'     regions of the CMDP matrix that are deemed negligible by an HMM.
#'     Default: `TRUE`
#' @param tau (`numeric` scalar) Tail loss probability used during HMM band
#'     calculation. Default: `1E-7`
#' @param mxsize (`numeric` scalar) The maximum DP matrix size, in Mb.
#'     Maximum potential memory usage is approximately cpu*mxsize, although
#'     this is usually not realized.
#'     See [Infernal](http://eddylab.org/infernal/) documentation for
#'     more information.
#' @param fixedtau (`logical` scalar) Turn off HMM band tightening.
#' @param maxtau (`numeric` scalar) Maximum allowed value for tau
#'     during band tightening. Default: `0.05`
#' @param small (`logical` scalar) Use the divide and conquer CYK alignment
#'     algorithm, greatly reducing memory consumption.
#' @param sfile (`character` filename) Dump per-sequence alignment score
#'     and timing information to file.
#' @param tfile (`character` filename) Dump tabular sequence tracebacks for
#'     each individual sequence to a file.
#' @param ifile (`character` filename) Dump per-sequence insert information
#'     to file.
#' @param elfile (`character` filename) Dump per-sequence EL state (local
#'     end) insert information to file.
#' @param mapali (`character` filename) Read the alignment from the file
#'     used
#'     to build the model aligns it as a single object to the CM, along with
#'     sequences in `"seq"`.
#' @param mapstr (`logical` scalar) Must be used in combination with
#'     `mapali`. Propogate structural information for any pseudoknots
#'     that exist in `mapali` to the output alignment.
#' @param dnaout (`logical` scalar) Output the alignments as DNA sequence
#'     alignments, instead of RNA ones.
#' @param noprob (`logical` scalar) Do not annotate the output alignment
#'     with posterior probabilities.
#' @param matchonly (`logical` scalar) Only include match columns in the
#'     output alignment, do not include any insertions relativeto the
#'     consensus model.
#' @param ileaved (`logical` scalar) Output the alignment in interleaved
#'     Stockholm format of a fixed width that may be more convenient for
#'     examination.
#' @param regress (`character` filename) Save an additional copy of the
#'     output alignment with no author information to file.
#' @param verbose (`logical` scalar) Output additional information in the
#'     tabular scores output.
#' @param cpu (`integer` scalar) The number of cpus to use.
#' @param mpi (`logical` scalar) Run as an MPI parallel program.
#' @param extra (`character` vector) Additional advanced options to pass
#'     to `cmalign`.
#' @param glocal (`logical` scalar) (Deprecated) see "global".
#'
#' @return the aligned sequences, as returned by [read_stockholm_msa()].
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
        nonneg_float_opt(mxsize),
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
            if (!requireNamespace("ShortRead", quietly = TRUE)) {
                stop("'seq' is a ShortRead object but the ShortRead package",
                     "is unavailable")
            }
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
    alnpipe <- pipe(args, open = "rt")
    out <- read_stockholm_msa(alnpipe, type = if (dnaout) "DNA" else "RNA")
    if (close(alnpipe) != 0) stop("cmalign had nonzero exit status.")
    out
}

#' Construct covariance model(s) from structually annotated alignment(s)
#'
#' Calls the standalone program `cmbuild` from the
#' [Infernal](http://eddylab.org/infernal/) package, which must be
#' installed.  For more information about options, see the
#' [Infernal documentation](http://eddylab.org/infernal/Userguide.pdf).
#'
#' @param msafile (`character` filename) filename of MSA file to read.
#' @param cmfile_out (`character` filename) filename to write CM to.
#' @param name (`character` string) name of the CM (option `-N` to
#'     cmbuild). The default uses the
#'     name(s) given in the alignment file, or the name of the alignment file.
#' @param force (`logical` scalar) overwrite `cmfile_out` if it exists.
#' @param summary_file (`character` filename) filename to write summary
#'     output to.
#' @param reannotated_msa (`character` filename) filename to write a
#'     reannotated alignment(s) to (option `-O` to cmbuild)
#' @param consensus_method (one of `"fast"`, `"hand"`, or
#'     `"noss"`) method to define consensus columns in alignment. (options
#'     `--fast`, `--hand`, and `--noss` to cmbuild)
#' @param symfrac (`integer` scalar) residue fraction threshold necessary
#'     to define a consensus column.
#' @param rsearch (`character` filename) RIBOSUM matrix to use in
#'     parameterizing emission scores.
#' @param null (`character` filename) null model file.
#' @param prior (`character` filename) Dirichlet prior file.
#' @param weights (`character`; one of `"wpb"`, `"wgsc"`, `"wnone"`, `"wgiven"`,
#'     or `"wblosum"`) sequence weighting method (options `--wpb`, `--wgsc`,
#'     `--wnone`, `--wgiven`, and `--wblosum` to cmbuild).
#' @param wid (`numeric` scalar) percent identity for clustering when
#'     `weights = "blosum"`.
#' @param eff_num (`character`; one of `"eent"` or `"enone"`) entropy weighting
#'     strategy (options `--eent` or `--enone` to cmbuild).
#' @param ere (`numeric` scalar) target mean match state relative entropy.
#' @param eminseq (`numeric` scalar) minimum allowed effective sequence number.
#' @param emaxseq (`numeric` scalar) maximum allowed effective sequence number.
#' @param ehmmre (`numeric` scalar) target HMM mean match state relative
#'     entropy.
#' @param eset (`numeric` scalar) effective sequence number for entropy
#'     weighting.
#' @param p7ere (`numeric` scalar) target mean match state relative entropy
#'     for the filter p7 HMM.
#' @param p7ml (`logical` scalar) use a maximum liklihood p7 HMM built from
#'     the CM.
#' @param refine (`character` filename) if given, the alignment is iteratively
#'     refined, and the final alignment is written to this file.
#' @param local (`logical` scalar) use local alignment with `refine` (option
#'     `-l` to cmbuild).
#' @param gibbs (`logical` scalar) use Gibbs sampling instead of maximum
#'     likelihood with `refine`.
#' @param seed (`integer` scalar) pseudorandom number seed for Gibbs sampling.
#' @param cyk (`logical` scalar) use the CYK alignment algorithm with `refine`.
#' @param notrunc (`logical` scalar) turn off the truncated alignment
#'     algorithm with `refine`.
#' @param extra (`character`) additional arguments to pass to cmbuild.
#' @param quiet (`logical` scalar) do not print cmbuild output to the console.
#'
#' @return `NULL`, invisibly
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

#' Calibrate a covariance model for searching
#'
#' @param cmfile (`character` file name) Filename of the CM to calibrate.
#' calibration data will be added to this file.
#' @param length (`numeric` scalar) Total length of random sequences to search,
#' in megabases. (parameter "`-L`" to `cmcalibrate`)
#' @param forecast (`logical` flag) Predict the running time of
#' the calibration of cmfile (with provided options) on the current machine and
#' exit. The calibration is not performed.
#' @param nforecast (`integer` count) With `forecast`, specify the number of
#' processors that will be used for the calibration.
#' @param memreq (`logical` flag) Predict the amount of required memory for
#' calibrating cmfile (with provided options) on the current machine and exit.
#' The calibration is not performed.
#' @param gtailn (`numeric` scalar) Fit the exponential tail for glocal Inside 
#' and glocal CYK to the `n` highest scores in the histogram tail, where `n` is
#' `gtailn` times `length` (the number of Mb searched).
#' @param ltailn (`numeric` scalar) Fit the exponential tail for local Inside 
#' and local CYK to the `n` highest scores in the histogram tail, where `n` is
#' `ltailn` times `length` (the number of Mb searched).
#' @param tailp (`numeric` scalar) Ignore `gtailn` and `ltailn` options and fit
#' the `tailp` fraction tail of the histogram to an exponential tail, for all
#' search modes.
#' @param hfile (`character` file name) File to save the histograms fit.
#' @param sfile (`character` file name) File to save the survival plot
#' information.
#' @param qqfile (`character` file name) File to save the quantile-quantile plot
#' information. 
#' @param ffile (`character` file name) File to save the statistics of
#' different exponential fits.
#' @param xfile (`character` file name) File to save a list of the scores in
#' each fit histogram.
#' @param seed (non-negative `integer` scalar) Random number seed. `0` seeds
#' the random number generator arbitrarily and non-reproducibly.
#' @param beta (non-negative `numeric` scalar) Beta tail loss probability used
#' for query-dependent banding (QDB).
#' @param nonbanded (`logical` flag) If `TRUE`, turn off QDB during E-value
#' calibration. This will slow down calibration.
#' @param nonull3 (`logical` flag) If `TRUE`, turn off the null3 post hoc
#' additional null model. This is not recommended unless you plan on using the
#' same option to `cmsearch` and/or `cmscan`.
#' @param random (`logical` flag) If `TRUE`, use the background null model of
#' the CM to generate the random sequences, instead of the more realistic HMM.
#' @param gc (`character` file name) A file to use to determine the nucleotide
#' frequencies for the random sequences.
#' @param cpu (`integer` count) The number of parallel workers to use.
#' @param mpi (`logical` flag) If `TRUE`, run as an MPI parallel program.
#' @param quiet (`logical` flag) If `TRUE`, suppress standard output of
#' `cmcalibrate`.
#'
#' @return NULL, invisibly
#' @export
cmcalibrate <- function(
    cmfile,
    length = 1.6,
    forecast = FALSE,
    nforecast = NULL,
    memreq = FALSE,
    gtailn = 250,
    ltailn = 750,
    tailp = NULL,
    hfile = NULL,
    sfile = NULL,
    qqfile = NULL,
    ffile = NULL,
    xfile = NULL,
    seed = 181L,
    beta = 1e-15,
    nonbanded = FALSE,
    nonull3 = FALSE,
    random = FALSE,
    gc = NULL,
    cpu = NULL,
    mpi = FALSE,
    quiet = FALSE
) {
    assertthat::assert_that(
        assertthat::is.readable(cmfile)
    )
    args <- c(
        nonneg_float_opt(length, "L"),
        flag_opt(forecast),
        count_opt(nforecast),
        flag_opt(memreq),
        nonneg_float_opt(gtailn),
        nonneg_float_opt(ltailn),
        nonneg_float_opt(tailp),
        string_opt(hfile),
        string_opt(sfile),
        string_opt(qqfile),
        string_opt(ffile),
        string_opt(xfile),
        count_opt(seed),
        nonneg_float_opt(beta),
        flag_opt(nonbanded),
        flag_opt(nonull3),
        flag_opt(random),
        infile_opt(gc),
        count_opt(cpu),
        flag_opt(mpi),
        cmfile
    )
    assertthat::assert_that(assertthat::is.flag(quiet))
    
    return <- system2(
        command = "cmcalibrate",
        args = args,
        stdout = if (isTRUE(quiet)) FALSE else "",
        stderr = if (isTRUE(quiet)) FALSE else ""
    )
    
    if (return != 0) stop("cmcalibrate failed with exit code ", return)
    invisible(NULL)
}