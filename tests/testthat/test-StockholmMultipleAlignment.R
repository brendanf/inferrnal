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
rna <- readRDS("msa_alignment.RDS")
rnama <-Biostrings::RNAMultipleAlignment(rna)

test_that("alternate initializations for RNAMultipleAlignment are identical", {
    expect_stockholm_equal(
        StockholmRNAMultipleAlignment(rna),
        StockholmRNAMultipleAlignment(rnama)
    )
    expect_stockholm_equal(
        StockholmRNAMultipleAlignment(rna),
        StockholmRNAMultipleAlignment(Biostrings::RNAStringSet(rna))
    )
})
test_that("initialization from RNAMultipleAlignment warns for extra arguments", {
    expect_warning(StockholmRNAMultipleAlignment(rnama, start = 0L))
    expect_warning(StockholmRNAMultipleAlignment(rnama, end = 0L))
    expect_warning(StockholmRNAMultipleAlignment(rnama, width = 0L))
    expect_warning(StockholmRNAMultipleAlignment(rnama, use.names = FALSE))
    expect_warning(StockholmRNAMultipleAlignment(rnama, rowmask = 3L))
    expect_warning(StockholmRNAMultipleAlignment(rnama, colmask = 3L))
})

dna <- chartr("U", "T", rna)
dnama <-Biostrings::DNAMultipleAlignment(dna)

test_that("alternate initializations for DNAMultipleAlignment are identical", {
    expect_stockholm_equal(
        StockholmDNAMultipleAlignment(dna),
        StockholmDNAMultipleAlignment(dnama)
    )
    expect_stockholm_equal(
        StockholmDNAMultipleAlignment(dna),
        StockholmDNAMultipleAlignment(Biostrings::DNAStringSet(dna))
    )
})
test_that("initialization from DNAMultipleAlignment warns for extra arguments", {
    expect_warning(StockholmDNAMultipleAlignment(dnama, start = 0L))
    expect_warning(StockholmDNAMultipleAlignment(dnama, end = 0L))
    expect_warning(StockholmDNAMultipleAlignment(dnama, width = 0L))
    expect_warning(StockholmDNAMultipleAlignment(dnama, use.names = FALSE))
    expect_warning(StockholmDNAMultipleAlignment(dnama, rowmask = 3L))
    expect_warning(StockholmDNAMultipleAlignment(dnama, colmask = 3L))
})

aa <- c(
    "O83071/192-246" = "MTCRAQLIAVPRASSLAE..AIACAQKM....RVSRVPVYERS",
    "O83071/259-312" = "MQHVSAPVFVFECTRLAY..VQHKLRAH....SRAVAIVLDEY",
    "O31698/18-71"   = "MIEADKVAHVQVGNNLEH..ALLVLTKT....GYTAIPVLDPS",
    "O31698/88-139"  = "EVMLTDIPRLHINDPIMK..GFGMVINN......GFVCVENDE",
    "O31699/88-139"  = "EVMLTDIPRLHINDPIMK..GFGMVINN......GFVCVENDE"
)
aama <- Biostrings::AAMultipleAlignment(aa)
test_that("initialization from AAMultipleAlignment warns for extra arguments", {
    expect_warning(StockholmAAMultipleAlignment(aama, start = 0L))
    expect_warning(StockholmAAMultipleAlignment(aama, end = 0L))
    expect_warning(StockholmAAMultipleAlignment(aama, width = 0L))
    expect_warning(StockholmAAMultipleAlignment(aama, use.names = FALSE))
    expect_warning(StockholmAAMultipleAlignment(aama, rowmask = 3L))
    expect_warning(StockholmAAMultipleAlignment(aama, colmask = 3L))
})

test_that("alternate initializations for AAMultipleAlignment are identical", {
    expect_stockholm_equal(
        StockholmAAMultipleAlignment(aa),
        StockholmAAMultipleAlignment(aama)
    )
    expect_stockholm_equal(
        StockholmAAMultipleAlignment(aa),
        StockholmAAMultipleAlignment(Biostrings::AAStringSet(aa))
    )
})
