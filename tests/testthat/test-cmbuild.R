test_cm <- tempfile("test", fileext = ".cm")

test_that("cmbuild can replicate RDP CM", {
    expect_null(
        cmbuild(
            msafile = system.file(
                file.path("extdata", "fungi_32S_LR5.stk"), package = "LSUx"
            ),
            cmfile_out = test_cm,
            force = TRUE,
            consensus_method = "hand",
            ere = 0.85,
            extra = "--iflank",
            verbose = TRUE)
    )
    # these lines have timestamps and system paths that are not expected to
    # replicate
  expect_identical(
      readLines(
          system.file(file.path("extdata", "fungi_32S_LR5.cm"), package = "LSUx")
      )[-c(2, 11, 12, 4328, 4337, 4338)],
      readLines(test_cm)[-c(2, 11, 12, 4328, 4337, 4338)]
  )
})
