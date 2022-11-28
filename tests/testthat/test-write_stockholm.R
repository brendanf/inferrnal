expect_stockholm_equal <- function(msa1, msa2) {
    expect_mapequal(as.character(msa1), as.character(msa2))
    if (length(msa1@GF) > 0 || length(msa2@GF) > 0)
        expect_mapequal(as.character(msa1@GF), as.character(msa2@GF))
    if (length(msa1@GC) > 0 || length(msa2@GC) > 0)
        expect_mapequal(as.character(msa1@GC), as.character(msa2@GC))
    expect_setequal(as.character(names(msa1@GR)), as.character(names(msa2@GR)))
    for (gr in names(msa1@GR)) {
        expect_mapequal(as.character(msa1@GR[[gr]]), as.character(msa2@GR[[gr]]))
    }
    expect_setequal(as.character(names(msa1@GS)), as.character(names(msa2@GS)))
    for (gs in names(msa1@GS)) {
        expect_mapequal(as.character(msa1@GS[[gs]]), as.character(msa2@GS[[gs]]))
    }
}

test_that("no change on write/read for AA", {
    t <- read_stockholm_msa("test.stk", type = "AA")
    f <- tempfile(fileext = ".stockholm")
    on.exit(unlink(f))
    writeStockholmMultipleAlignment(t, f)
    t2 <- read_stockholm_msa(f, "AA")
    expect_stockholm_equal(t, t2)
})

test_that("no change on write/read for RNA", {
    t <- read_stockholm_msa(system.file("extdata/sample.stk", package="inferrnal"), type = "RNA")
    f <- tempfile(fileext = ".stockholm")
    on.exit(unlink(f))
    writeStockholmMultipleAlignment(t, f)
    t2 <- read_stockholm_msa(f, "RNA")
    expect_stockholm_equal(t, t2)
})
