check_stockholm_msa <- function(object) {
    errors <- character()
    width <- ncol(object)
    for (gs in names(object@GS)) {
        mismatch <- setdiff(names(object@GS[[gs]]), names(object@unmasked))
        if (length(mismatch)) {
            msg <- paste("sequence", mismatch,
                         "found in GS", gs, "annotation but not alignment.")
            errors <- c(errors, msg)
        }
    }
    for (gr in names(object@GR)) {
        mismatch <- setdiff(names(object@GR[[gr]]), names(object@unmasked))
        if (length(mismatch) > 0) {
            msg <- paste("sequence", mismatch,
                         "found in GR", gr, "annotation but not alignment.")
            errors <- c(errors, msg)
        }
        wronglength <- which(Biostrings::nchar(object@GR[[gr]]) != width)
        if (length(wronglength) > 0) {
            msg <- paste("GR", gr, "annotation for sequence",
                         names(object@GR[[gr]])[wronglength], "has",
                         Biostrings::nchar(object@GR[[gr]])[wronglength],
                         "characters but alignment width is", width)
            errors <- c(errors, msg)
        }
    }
    
    wronglength <- which(Biostrings::nchar(object@GC) != width)
    if (length(wronglength) > 0) {
        msg <- paste("GC", names(object@GC)[wronglength], "annotation has",
                     Biostrings::nchar(object@GC)[wronglength],
                     "characters but alignment width is", width)
        errors <- c(errors, msg)
    }
    if (length(errors) == 0) TRUE else errors
}


#' StockholmMultipleAlignment objects
#' 
#' The `StockholmMultipleAlignment` class contains a multiple sequence
#' alignment along with its annotations, as defined for the Stockholm file
#' format.
#' 
#' Although the `StockholmMultipleAlignment` class is agnostic about the
#' specific tags used, the following tags are the most likely to be recognized
#' by Infernal or other software which reads or writes Stockholm files:
#' 
#' \tabular{lll}{
#' Type \tab Tag \tab Description \cr
#' GF \tab ID \tab IDentifier \cr
#' GF \tab AC \tab ACcession \cr
#' GF \tab DE \tab DEscription \cr
#' GF \tab AU \tab AUthor \cr
#' GF \tab GA \tab GAthering threshold \cr
#' GF \tab NC \tab Noise Cutoff \cr
#' GF \tab TC \tab Trusted Cutoff \cr
#' GS \tab WT \tab WeighT \cr
#' GS \tab AC \tab ACcession number \cr
#' GS \tab DE \tab DEscription \cr
#' GS \tab DR \tab Database Reference \cr
#' GS \tab OS \tab OrganiSm (species) \cr
#' GS \tab OC \tab Organism Classification (clade, etc.) \cr
#' GS \tab LO \tab Look (Color, etc.) \cr
#' GR \tab SS \tab Secondary Structure \cr
#' GR \tab SA \tab Surface Accessibility \cr
#' GR \tab TM \tab TransMembrane \cr
#' GR \tab PP \tab Posterior Probability \cr
#' GR \tab LI \tab LIgand binding \cr
#' GR \tab AS \tab Active Site \cr
#' GR \tab pAS \tab AS - Pfam predicted \cr
#' GR \tab sAS \tab AS - from SwissProt \cr
#' GR \tab IN \tab INtron (in or after) \cr
#' GC \tab RF \tab ReFerence \cr
#' GC \tab SS_cons \tab Secondary Structure consensus \cr
#' GC \tab SA_cons \tab Surface Accessibility consensus \cr
#' GC \tab TM_cons \tab TransMembrane consensus \cr
#' GC \tab PP_cons \tab Posterior Probability consensus \cr
#' GC \tab LI_cons \tab LIgand binding consensus \cr
#' GC \tab AS_cons \tab Active Site consensus \cr
#' GC \tab pAS_cons \tab AS - Pfam predicted consensus \cr
#' GC \tab sAS_cons \tab AS - from SwissProt consensus \cr
#' GC \tab IN_cons \tab INtron (in or after) consensus \cr
#' }
#'
#' @slot GF [`BStringSet`][Biostrings::XStringSet-class]. Free-text annotations
#' which belong to the alignment file as a whole.  The name of each element is
#' a tag identifying the type of data. (See Details).
#' @slot GS [`BStringSetList`][Biostrings::XStringSetList-class]. Free-text
#' annotations which belong to the individual sequences in the alignment. The
#' name of each [`BStringSet`][Biostrings::XStringSet-class] is a tag
#' identifying the type of data. (See Details). Names of individual
#' [`BString`][Biostrings::XString-class] elements match the names of sequences
#' in the alignment, but there is no requirement that every sequence must be
#' annotated for every tag.
#' @slot GR [`BStringSetList`][Biostrings::XStringSetList-class]. Annotations
#' for individual residues in the alignment. The name of each
#' [`BStringSet`][Biostrings::XStringSet-class] is a tag identifying the type
#' of data. (See Details). Names of individual
#' [`BString`][Biostrings::XString-class] elements match the names of sequences
#' in the alignment, but there is no requirement that every sequence must be
#' annotated for every tag. Unlike GS tags, the width of all elements must be
#' the same, and must match the width of the alignment.
#' @slot GC [`BStringSet`][Biostrings::XStringSet-class]. Annotations
#' which belong to each column of the alignment as a whole. The name of each
#' element is a tag identifying the type of data. (See Details).  Unlike GF
#' tags, the width of all elements must be the same, and must match the width
#' of the alignment.
#'
#' @export
setClass(
    Class = "StockholmMultipleAlignment",
    slots = list(
        GF = "BStringSet",
        GS = "BStringSetList",
        GR = "BStringSetList",
        GC = "BStringSet"
    ),
    validity = check_stockholm_msa,
    prototype = list(
        GF = Biostrings::BStringSet(),
        GS = Biostrings::BStringSetList(),
        GR = Biostrings::BStringSetList(),
        GC = Biostrings::BStringSet()
    ),
    contains = c("VIRTUAL")
)

setClass(
    Class = "StockholmDNAMultipleAlignment",
    contains = c("StockholmMultipleAlignment", "DNAMultipleAlignment")
) -> .StockholmDNAMultipleAlignment

setClass(
    Class = "StockholmRNAMultipleAlignment",
    contains = c("StockholmMultipleAlignment", "RNAMultipleAlignment")
) -> .StockholmRNAMultipleAlignment

setClass(
    Class = "StockholmAAMultipleAlignment",
    contains = c("StockholmMultipleAlignment", "AAMultipleAlignment")
) -> .StockholmAAMultipleAlignment

setMethod(
    initialize,
    "StockholmMultipleAlignment",
    function(.Object, GF = character(), GS = list(),
             GR = list(), GC = character(), ...) {
        callNextMethod(.Object, ..., GF = GF, GS = GS, GR = GR, GC = GC)
    }
)

nonmissing_warning <- function(arg, type = c("RNA", "DNA", "AA")) {
    type = match.arg(type)
    warning("'", arg, "' is ignored when `x` is not `character` or `",
            type, "StringSet`")
}

#' @rdname StockholmMultipleAlignment-class
#' @param x (`character` vector, aligned
#' [`XStringSet`][Biostrings::XStringSet-class] or
#' [`MultipleAlignment`][Biostrings::MultipleAlignment-class] of an appropriate
#' type) multiple sequence alignment without Stockholm extensions
#' @param start,end,width,use.names,rowmask,colmask passed to the appropriate
#' [`MultipleAlignment`][Biostrings::MultipleAlignment-class] constructor,
#' unless "x" is already a
#' [`MultipleAlignment`][Biostrings::MultipleAlignment-class]
#' @param GF (named `character` vector or
#' [`BStringSet`][Biostrings::XStringSet-class]) Free-text annotations
#' which belong to the alignment file as a whole.  The name of each element is
#' a tag identifying the type of data. (See Details).
#' @param GS (named `list` of named `character` vectors or
#' [`BStringSetList`][Biostrings::XStringSetList-class]) Free-text
#' annotations which belong to the individual sequences in the alignment. The
#' names of the outer `list` or
#' [`BStringSetList`][Biostrings::XStringSetList-class] are tags
#' identifying the type of data for each element. (See Details). Names of inner
#' `character` vectors or [`BStringSet`][Biostrings::XStringSet-class]s match
#' the names of sequences in the alignment, but there is no requirement that
#' every sequence must be annotated for every tag.
#' @param GR (named `list` of named `character` vectors or
#' [`BStringSetList`][Biostrings::XStringSetList-class]) Annotations
#' for individual residues in the alignment. The names of the outer `list` or
#' [`BStringSetList`][Biostrings::XStringSetList-class] are tags
#' identifying the type of data for each element. (See Details). Names of inner
#' `character` vectors or [`BStringSet`][Biostrings::XStringSet-class]s match
#' the names of sequences in the alignment, but there is no requirement that
#' every sequence must be annotated for every tag. Unlike GS tags, the width of
#' all elements must be the same, and must match the width of the alignment.
#' @param GC (named `character` vector or
#' [`BStringSet`][Biostrings::XStringSet-class]) Annotations
#' which belong to each column of the alignment as a whole. The name of each
#' element is a tag identifying the type of data. (See Details).  Unlike GF
#' tags, the width of all elements must be the same, and must match the width
#' of the alignment.
#' @export
#' @return a new `StockholmMultipleAlignment` object
#' 
#' @examples
#' # Typically a StockholmMultipleAlignment object is read from a file created
#' # by other software, but it can also be created manually.
#' # This example reproduces the example file given in the Stockholm format
#' # definition.
#' samp <- StockholmAAMultipleAlignment(
#'     x = c(
#'         "O83071/192-246" = "MTCRAQLIAVPRASSLAE..AIACAQKM....RVSRVPVYERS",
#'         "O83071/259-312" = "MQHVSAPVFVFECTRLAY..VQHKLRAH....SRAVAIVLDEY",
#'         "O31698/18-71"   = "MIEADKVAHVQVGNNLEH..ALLVLTKT....GYTAIPVLDPS",
#'         "O31698/88-139"  = "EVMLTDIPRLHINDPIMK..GFGMVINN......GFVCVENDE",
#'         "O31699/88-139"  = "EVMLTDIPRLHINDPIMK..GFGMVINN......GFVCVENDE"
#'     ),
#'     GF = c(
#'         ID = "CBS",
#'         AC = "PF00571",
#'         AU = "Bateman A",
#'         CC = paste("CBS domains are small intracellular modules mostly",
#'                    "found in 2 or four copies within a protein."),
#'         SQ = "67"
#'     ),
#'     GS = list(
#'         # ACcession number
#'         AC = c(
#'             "O31698/18-71" = "O31698",
#'             "O83071/192-246" = "O83071",
#'             "O83071/259-312" = "O83071",
#'             "O31698/88-139" = "O31698"
#'         ),
#'         # OrganiSm
#'         OS = c("O31698/88-139" = "Bacillus subtilis")
#'     ),
#'     GR = list(
#'         # Surface Accessibility
#'         SA = c(
#'             "O83071/192-246" = "999887756453524252..55152525....36463774777"
#'         ),
#'         # Secondary Structure
#'         SS = c(
#'             "O83071/259-312" = "CCCCCHHHHHHHHHHHHH..EEEEEEEE....EEEEEEEEEEE",
#'             "O31698/18-71"   = "CCCHHHHHHHHHHHHHHH..EEEEEEEE....EEEEEEEEHHH",
#'             "O31698/88-139"  = "CCCCCCCHHHHHHHHHHH..HEEEEEEE....EEEEEEEEEEH"
#'         ),
#'         # Active Site
#'         AS = c(
#'             "O31699/88-139"  = "________________*__________________________"
#'         ),
#'         # INtron
#'         IN = c(
#'             "O31699/88-139"  = "____________1______________2__________0____"
#'         )
#'     ),
#'     GC = c(
#'         # Secondary Structure consensus
#'         SS_cons = "CCCCCHHHHHHHHHHHHH..EEEEEEEE....EEEEEEEEEEH"
#'     )
#'  )
#'  samp
StockholmDNAMultipleAlignment <- function(
    x = character(),
    start = NA,
    end = NA,
    width = NA,
    use.names = TRUE,
    rowmask = NULL,
    colmask = NULL,
    GF = character(),
    GS = list(),
    GR = list(),
    GC = character()
){
    if (methods::is(x, "character") || methods::is(x, "DNAStringSet")) {
        x <- Biostrings::DNAMultipleAlignment(
            x = x,
            start = start,
            end = end,
            width = width,
            use.names = use.names,
            rowmask = rowmask,
            colmask = colmask
        )
    } else if (methods::is(x, "DNAMultipleAlignment")) {
        if (!all(is.na(start))) nonmissing_warning("start", "DNA")
        if (!all(is.na(end))) nonmissing_warning("end", "DNA")
        if (!all(is.na(width))) nonmissing_warning("width", "DNA")
        if (!isTRUE(use.names)) nonmissing_warning("use.names", "DNA")
        if (!is.null(rowmask)) nonmissing_warning("rowmask", "DNA")
        if (!is.null(colmask)) nonmissing_warning("colmask", "DNA")
    } else {
        stop(
            "'x' should be `character`, `DNAStringSet`, ",
            "or `DNAMultipleAlignment`"
        )
    }
    if (methods::is(GF, "character")) GF <- Biostrings::BStringSet(GF)
    if (methods::is(GS, "list")) GS <-
            Biostrings::BStringSetList(lapply(GS, Biostrings::BStringSet))
    if (methods::is(GR, "list")) GR <-
            Biostrings::BStringSetList(lapply(GR, Biostrings::BStringSet))
    if (methods::is(GC, "character")) GC <- Biostrings::BStringSet(GC)
    methods::new(
        "StockholmDNAMultipleAlignment",
        x,
        GF = GF,
        GS = GS,
        GR = GR,
        GC = GC
    )
}

#' @rdname StockholmMultipleAlignment-class
#' @export
StockholmRNAMultipleAlignment <- function(
    x = character(),
    start = NA,
    end = NA,
    width = NA,
    use.names = TRUE,
    rowmask = NULL,
    colmask = NULL,
    GF = character(),
    GS = list(),
    GR = list(),
    GC = character()
){
    if (methods::is(x, "character") || methods::is(x, "RNAStringSet")) {
        x <- Biostrings::RNAMultipleAlignment(
            x = x,
            start = start,
            end = end,
            width = width,
            use.names = use.names,
            rowmask = rowmask,
            colmask = colmask
        )
    } else if (methods::is(x, "RNAMultipleAlignment")) {
        if (!all(is.na(start))) nonmissing_warning("start", "RNA")
        if (!all(is.na(end))) nonmissing_warning("end", "RNA")
        if (!all(is.na(width))) nonmissing_warning("width", "RNA")
        if (!isTRUE(use.names)) nonmissing_warning("use.names", "RNA")
        if (!is.null(rowmask)) nonmissing_warning("rowmask", "RNA")
        if (!is.null(colmask)) nonmissing_warning("colmask", "RNA")
    } else {
        stop(
            "'x' should be `character`, `RNAStringSet`, ",
            "or `RNAMultipleAlignment`"
        )
    }
    if (methods::is(GF, "character")) GF <- Biostrings::BStringSet(GF)
    if (methods::is(GS, "list")) GS <-
            Biostrings::BStringSetList(lapply(GS, Biostrings::BStringSet))
    if (methods::is(GR, "list")) GR <-
            Biostrings::BStringSetList(lapply(GR, Biostrings::BStringSet))
    if (methods::is(GC, "character")) GC <- Biostrings::BStringSet(GC)
    methods::new(
        "StockholmRNAMultipleAlignment",
        x,
        GF = GF,
        GS = GS,
        GR = GR,
        GC = GC
    )
}

#' @rdname StockholmMultipleAlignment-class
#' @export
StockholmAAMultipleAlignment <- function(
    x = character(),
    start = NA,
    end = NA,
    width = NA,
    use.names = TRUE,
    rowmask = NULL,
    colmask = NULL,
    GF = character(),
    GS = list(),
    GR = list(),
    GC = character()
){
    if (methods::is(x, "character") || methods::is(x, "AAStringSet")) {
        x <- Biostrings::AAMultipleAlignment(
            x = x,
            start = start,
            end = end,
            width = width,
            use.names = use.names,
            rowmask = rowmask,
            colmask = colmask
        )
    } else if (methods::is(x, "AAMultipleAlignment")) {
        if (!all(is.na(start))) nonmissing_warning("start", "AA")
        if (!all(is.na(end))) nonmissing_warning("end", "AA")
        if (!all(is.na(width))) nonmissing_warning("width", "AA")
        if (!isTRUE(use.names)) nonmissing_warning("use.names", "AA")
        if (!is.null(rowmask)) nonmissing_warning("rowmask", "AA")
        if (!is.null(colmask)) nonmissing_warning("colmask", "AA")
    } else {
        stop(
            "'x' should be `character`, `AAStringSet`, ",
            "or `AAMultipleAlignment`"
        )
    }
    if (methods::is(GF, "character")) GF <- Biostrings::BStringSet(GF)
    if (methods::is(GS, "list")) GS <-
            Biostrings::BStringSetList(lapply(GS, Biostrings::BStringSet))
    if (methods::is(GR, "list")) GR <-
            Biostrings::BStringSetList(lapply(GR, Biostrings::BStringSet))
    if (methods::is(GC, "character")) GC <- Biostrings::BStringSet(GC)
    methods::new(
        "StockholmAAMultipleAlignment",
        x,
        GF = GF,
        GS = GS,
        GR = GR,
        GC = GC
    )
}

show_stockholm_msa <- function(object) {
    callNextMethod(object)
    if (length(object@GF) > 0) {
        cat("\nGF (file) annotations:\n")
        show(object@GF)
    }
    if (length(object@GS) > 0) {
        for (gs in names(object@GS)) {
            cat("\nGS (sequence)", gs, "annotations:\n")
            show(object@GS[[gs]])
        }
    }
    if (length(object@GR) > 0) {
        for (gr in names(object@GR)) {
            cat("\nGR (residue)", gr, "annotations:\n")
            show(object@GR[[gr]])
        }
    }
    if (length(object@GC) > 0) {
        cat("\nGC (column) annotations:\n")
        show(object@GC)
    }
}

#' @importFrom methods callNextMethod new
#' @importMethodsFrom methods show
setMethod(show, "StockholmDNAMultipleAlignment", show_stockholm_msa)
setMethod(show, "StockholmRNAMultipleAlignment", show_stockholm_msa)
setMethod(show, "StockholmAAMultipleAlignment", show_stockholm_msa)

narrow_stockholm_msa <- function(x, start = NA, end = NA, width = NA, use.names = TRUE) {
    new <- x
    sew <- IRanges::solveUserSEW(ncol(new), start = start, end = end, width = width)
    new@colmask <- IRanges::restrict(new@colmask, start = sew@start, end = sew@start + sew@width - 1L, use.names = use.names)
    new@colmask <- IRanges::shift(new@colmask, shift = 1L - sew@start, use.names = use.names)
    new@unmasked <- IRanges::narrow(new@unmasked, start = sew@start, width = sew@width, use.names = use.names)
    new@GC <- IRanges::narrow(new@GC, start = sew@start, width = sew@width, use.names = use.names)
    for (gr in names(new@GR)) {
        new@GR[[gr]] <- IRanges::narrow(new@GR[[gr]], start = sew@start, width = sew@width, use.names = TRUE)
    }
    methods::validObject(new)
    new
}

#' @importMethodsFrom IRanges narrow
setMethod(narrow, signature(x = "StockholmMultipleAlignment"), narrow_stockholm_msa)