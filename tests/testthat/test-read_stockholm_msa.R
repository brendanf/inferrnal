test_that("fails on integer", {
  expect_error(read_stockholm_msa(1L))
})

test_that("fails on double", {
    expect_error(read_stockholm_msa(1))
})

test_that("fails on logical", {
    expect_error(read_stockholm_msa(FALSE))
})

test_that("fails on nonfile", {
    expect_error(read_stockholm_msa(tempfile("fake")))
})

ref <- list(
    alignment = Biostrings::RNAMultipleAlignment(readRDS("msa_alignment.RDS")),
    GF = readRDS("msa_GF.RDS"),
    GS = list(),
    GR = list(
        PP = Biostrings::BStringSet(readRDS("msa_GR_PP.RDS"))
    ),
    GC = list(
        SS_cons = Biostrings::BString(readRDS("msa_GC_SS_cons.RDS")),
        RF = Biostrings::BString(readRDS("msa_GC_RF.RDS"))
    )
)


test_that("can read file", {
    expect_equal(
        read_stockholm_msa(
            system.file(
                file.path("extdata", "sample.stk"),
                package = "inferrnal"
            )
        ),
        ref
    )
})

test_that("reads example from format definition", {
  expect_snapshot_output(read_stockholm_msa("test.stk", type = "aa"))
})
