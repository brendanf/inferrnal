#' Download a covariance model from RFAM
#'
#' @param cm (\code{character} string) the name of the CM to download. Typically
#'     this is in the format \code{"RFxxxxx"} where the \code{"x"}'s represent
#'     decimal digits.
#' @param filename (\code{character} string) file to save the CM in.  The
#'     default is to save the CM in the current working directory in a file
#'     named \code{"RFxxxxx.cm"}.
#' @param ... additional arguments passed on to
#'     \code{\link[utils]{download.file}}.
#'
#' @return The name of the downloaded file, invisibly.
#' @export
#'
download_rfam_cm <- function(cm, filename = paste0(cm, ".cm"), ...) {
    assertthat::assert_that(
        assertthat::is.string(cm),
        assertthat::is.string(filename)
    )
    utils::download.file(
        paste0("https://rfam.xfam.org/family/", cm, "/cm"),
        destfile = filename,
        ...
    )
    invisible(filename)
}
