stockholm_test <- c(
    "testtaxon AGCG--AGC",
    "#=GC test -(((..)))"
)


test_that("parse_stockholm_msa_chunk works for a sequence", {
    expect_equal(
        inferrnal:::parse_stockholm_msa_chunk(stockholm_test[1], pos = 1, acc = list()),
        list(testtaxon = "AGCG--AGC")
    )
})

test_that("parse_stockholm_msa_chunk works for a GC line", {
    expect_equal(
        inferrnal:::parse_stockholm_msa_chunk(stockholm_test[2], pos = 1, acc = list()),
        structure(list(), test = "-(((..)))")
    )
})

test_that("parse_stockholm_msa_chunk works for a sequence + GC line", {
    expect_equal(
        inferrnal:::parse_stockholm_msa_chunk(stockholm_test, pos = 1, acc = list()),
        structure(list(testtaxon = "AGCG--AGC"), test = "-(((..)))")
    )
})
