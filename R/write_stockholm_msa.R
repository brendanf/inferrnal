#' Write an annotated alignment in Stockholm format
#'
#' @param x ([`StockholmMultipleAlignment`][StockholmMultipleAlignment-class])
#'     the object to write to the file
#' @param connection ([`connection`][base::connections] object or file name)
#'     the file to write the alignment to
#'
#' @return the connection, invisibly
#' @export
#'
#' @examples
#' # read an example rRNA alignment
#' msafile <- sample_rRNA_stk()
#' msa <- read_stockholm_msa(msafile)
#' msa
#' 
#' # write to another file
#' temp <- tempfile(fileext=".stk")
#' writeStockholmMultipleAlignment(msa, temp)
writeStockholmMultipleAlignment <- function(x, connection) {
  assertthat::assert_that(
    methods::is(x, "StockholmDNAMultipleAlignment") ||
      methods::is(x, "StockholmRNAMultipleAlignment") ||
      methods::is(x, "StockholmAAMultipleAlignment")
  )
  if (methods::is(connection, "connection") && !isOpen(connection)) {
    open(connection, "w")
    on.exit(close(connection))
  }
  cat("# STOCKHOLM 1.0\n", file = connection)
  if (length(x@GF) > 0) {
      cat(
          paste0(
              "#=GF ", names(x@GF), " ", as.character(x@GF), "\n", collapse=""
          ),
          file = connection,
          append = TRUE
      )
  }
  if (length(x@GS) > 0) {
      for (gs in names(x@GS)) {
          cat(
              paste(
                  "#=GS", names(x@GS[[gs]]), gs, as.character(x@GS[[gs]]),
                  collapse = "\n"
              ),
              "\n",
              sep = "",
              file = connection,
              append = TRUE
          )
      }
  }
  seqnames <- names(x@unmasked)
  gr_names_max <- max(
      0,
      vapply(
          names(x@GR),
          function(tag) {
              max(nchar(names(x@GR[[tag]]))) + 6L + nchar(tag)
          },
          3L
      )
  )
  gc_names_max <- max(nchar(names(x@GC)), 0) + 5L
  max_w <- max(nchar(seqnames), gr_names_max, gc_names_max)
  for (s in seqnames) {
      cat(
          s, strrep(" ", max_w - nchar(s) + 1L),
          as.character(x@unmasked[s]), "\n",
          sep = "",
          file = connection,
          append = TRUE
      )
      for (gr in names(x@GR)) {
          if (s %in% names(x@GR[[gr]])) {
              cat(
                  "#=GR ", s, " ", gr,
                  strrep(" ", max_w - nchar(s) - nchar(gr) - 5L),
                  as.character(x@GR[[gr]][s]), "\n",
                  sep = "",
                  file = connection,
                  append = TRUE
              )
          }
      }
  }
  if (length(x@GC) > 0) {
      for (gc in names(x@GC)) {
          cat("#=GC ", gc, strrep(" ", max_w - nchar(gc) - 4L),
                  as.character(x@GC[[gc]]), "\n",
              sep = "",
              file = connection,
              append = TRUE
          )
      }
      cat("//\n", file = connection, append = TRUE)
  }
  invisible(connection)
}