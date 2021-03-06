test_cm <- tempfile("test", fileext = ".cm")

test_that("cmbuild can replicate RDP CM", {
    expect_null(
        cmbuild(
            msafile = stk_5_8S(),
            cmfile_out = test_cm,
            force = TRUE,
            quiet = TRUE)
    )
    # these lines have timestamps and system paths that are not expected to
    # replicate, or are for the calibrated cm.
    # also one line has a deviation in the least significant digit.
    infernal_ver <- gsub(
        "# INFERNAL ([0-9.]+).+",
        "\\1",
        system2("cmbuild", "-h", stdout = TRUE)[2]
    )
    skip_if(utils::compareVersion(infernal_ver,"1.1.3") >= 0,
            "skipping comparison with RFAM CM due to changes in Infernal 1.1.3")
    expect_identical(
        readLines(cm_5_8S())[-c(1, 13, 14, 15, 32:35, 662, 665, 674:676, 913)],
        readLines(test_cm)[-c(1, 13, 14, 657, 660, 669:670, 907)]
    )
})
