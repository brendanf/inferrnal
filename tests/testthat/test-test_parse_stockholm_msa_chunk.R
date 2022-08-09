stockholm_test <- c(
    "#=GF AU Fake Author",
    "testtaxon AGCG--AGC",
    "#=GR testtaxon PP *********",
    "#=GS testtaxon AC 000000",
    "#=GC test -(((..)))",
    "#=GF AU et al.",
    "#=GS testtaxon AC 000001"
)


test_that("parse_stockholm_msa_chunk works for a sequence", {
    expect_equal(
        inferrnal:::parse_stockholm_msa_chunk(stockholm_test[2], pos = 1, acc = list()),
        list(alignment = c(testtaxon = "AGCG--AGC"))
    )
})

test_that("parse_stockholm_msa_chunk works for a GC line", {
    expect_equal(
        inferrnal:::parse_stockholm_msa_chunk(stockholm_test[5], pos = 1, acc = list()),
        list(GC = c(test = "-(((..)))"))
    )
})

test_that("parse_stockholm_msa_chunk works for a GS line", {
    expect_equal(
        inferrnal:::parse_stockholm_msa_chunk(stockholm_test[4], pos = 1, acc = list()),
        list(GS = list(AC = c(testtaxon = "000000")))
    )
})

test_that("parse_stockholm_msa_chunk works for split GS line", {
    expect_equal(
        inferrnal:::parse_stockholm_msa_chunk(stockholm_test[c(4, 7)], pos = 1, acc = list()),
        list(GS = list(AC = c(testtaxon = "000000 000001")))
    )
})

test_that("parse_stockholm_msa_chunk works for a GF line", {
    expect_equal(
        inferrnal:::parse_stockholm_msa_chunk(stockholm_test[1], pos = 1, acc = list()),
        list(GF = c(AU = "Fake Author"))
    )
})

test_that("parse_stockholm_msa_chunk works with split GF line", {
    expect_equal(
        inferrnal:::parse_stockholm_msa_chunk(stockholm_test[c(1, 6)], pos = 1, acc = list()),
        list(GF = c(AU = "Fake Author et al."))
    )
})

test_that("parse_stockholm_msa_chunk works for a GR line", {
    expect_equal(
        inferrnal:::parse_stockholm_msa_chunk(stockholm_test[3], pos = 1, acc = list()),
        list(GR = list(PP = c(testtaxon = "*********")))
    )
})

test_that("parse_stockholm_msa_chunk works for a sequence with annotations", {
    expect_equal(
        inferrnal:::parse_stockholm_msa_chunk(stockholm_test, pos = 1, acc = list()),
        list(
            GF = c(AU = "Fake Author et al."),
            GC = c(test = "-(((..)))"),
            GS = list(
                AC = c(testtaxon = "000000 000001")
                ),
            GR = list(
                PP = c(testtaxon = "*********")
                ),
            alignment = c(testtaxon = "AGCG--AGC")
        )
    )
})
