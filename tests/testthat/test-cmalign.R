expect_stockholm_equal <- function(msa1, msa2) {
    expect_mapequal(as.character(msa1), as.character(msa2))
    expect_mapequal(as.character(msa1@GF), as.character(msa2@GF))
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

ref <- StockholmRNAMultipleAlignment(
    Biostrings::RNAMultipleAlignment(readRDS("cmalign_alignment.RDS")),
    GF = c(
        AU = gsub(
            "# INFERNAL (.+) \\(.+\\)",
            "Infernal \\1",
            system("cmalign -h", intern = TRUE)[2]
        )
    ),
    GS = list(),
    GR = list(
        PP = readRDS("cmalign_GR_PP.RDS")
    ),
    GC = c(
        SS_cons = readRDS("cmalign_GC_SS_cons.RDS"),
        RF = readRDS("cmalign_GC_RF.RDS")
    )
)

test_that("cmalign works without regression", {
    expect_stockholm_equal(ref, cmalign(cm_5_8S(), sample_rRNA_5_8S(), cpu = 1))
    
})

test_that("glocal is deprecated in cmalign", {
    expect_warning(cmalign(cm_5_8S(), sample_rRNA_5_8S(), cpu = 1,
                           glocal = TRUE))
})

test_that("glocal and global don't work together", {
    expect_error(
        suppressWarnings(
            cmalign(cm_5_8S(), sample_rRNA_5_8S(), cpu = 1,
                    glocal = TRUE, global = TRUE)
        )
    )
})

rnass <- Biostrings::readRNAStringSet(sample_rRNA_5_8S())
dnass <- Biostrings::DNAStringSet(rnass)
rnachar <- as.character(rnass)
dnachar <- as.character(dnass)
dnasr <- ShortRead::ShortRead(dnass, Biostrings::BStringSet(names(dnass)))

test_that("cmalign works for different input formats", {
    expect_stockholm_equal(cmalign(cm_5_8S(), rnass, cpu = 1), ref)
    expect_stockholm_equal(cmalign(cm_5_8S(), dnass, cpu = 1), ref)
    expect_stockholm_equal(cmalign(cm_5_8S(), rnachar, cpu = 1), ref)
    expect_stockholm_equal(cmalign(cm_5_8S(), dnachar, cpu = 1), ref)
    expect_stockholm_equal(cmalign(cm_5_8S(), dnasr, cpu = 1), ref)
    expect_error(cmalign(cm_5_8S(), "This is a bogus sequence file", cpu = 1))
    expect_error(cmalign(cm_5_8S(), 17, cpu = 1))
})

test_that("error when cmalign fails", {
    expect_error(cmalign(cm_5_8S(), sample_rRNA_fasta(), cpu = 1, mxsize = 0.01))
})
