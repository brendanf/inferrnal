

write_stockholm_msa <- function(alignment, GF, GS, GR, GC, connection) {
  assertthat::assert_that(
    methods::is(alignment, "DNAMultipleAlignment") ||
      methods::is(alignment, "RNAMultipleAlignment") ||
      methods::is(alignment, "AAMultipleAlignment"),
    is.character(GF),
    is.list(GS),
    all(vapply(GS, methods::is, TRUE, "BStringSet")),
    is.list(GR),
    all(vapply(GR, methods::is, TRUE, "BStringSet")),
    is.list(GC),
    all(vapply(GC, methods::is, TRUE, "BString"))
  )
  if (methods::is(connection, "connection") && !isOpen(connection)) {
    open(connection, "w")
    on.exit(close(connection))
  }
  cat("# STOCKHOLM 1.0\n", file = connection)
  cat(
    paste0("#=GF ", names(GF), " ", GF, "\n", sep = " ", collapse=""),
    file = connection
  )
  for (gs in names(GS)) {
    cat(
      paste0("#=GS ", names(GS[[gs]]), " ", gs, " ", as.character(GS[[gs]]), "\n"),
      file = connection
    )
  }
  seqnames <- names(alignment@unmasked)
  gr_names_max <- max(
    vapply(
      names(GR),
      function(tag) {
        max(nchar(names(GR[[tag]]))) + 6L + nchar(tag)
      },
      3L
    )
  )
  gc_names_max <- max(nchar(names(GC))) + 5L
  max_w <- max(nchar(seqnames), gr_names_max, gc_names_max)
  for (s in seqnames) {
    cat(
      s, strrep(" ", max_w - nchar(s) + 1L),
      as.character(alignment@unmasked[s]), "\n",
      file = connection
    )
    for (gr in names(GR)) {
      cat(
        "#=GR ", s, " ", gr,
        strrep(" ", max_w - nchar(s) - nchar(gr) - 5L),
        as.character(GR[[gr]][s]), "\n",
        file = connection
      )
    }
  }
  for (gc in names(GC)) {
    cat("#=GC ", gc, strrep(" ", max_w - nchar(gc) - 4L),
        as.character(GC[[gc]]), "\n",
        file = connection
    )
  }
  cat("//\n", file = connection)
}