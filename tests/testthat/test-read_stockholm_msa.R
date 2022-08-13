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

ref <- StockholmRNAMultipleAlignment(
    Biostrings::RNAMultipleAlignment(readRDS("msa_alignment.RDS")),
    GF = readRDS("msa_GF.RDS"),
    GS = list(),
    GR = list(
        PP = readRDS("msa_GR_PP.RDS")
    ),
    GC = c(
        SS_cons = readRDS("msa_GC_SS_cons.RDS"),
        RF = readRDS("msa_GC_RF.RDS")
    )
)


test_that("can read file", {
    expect_stockholm_equal(
        read_stockholm_msa(
            system.file(
                file.path("extdata", "sample.stk"),
                package = "inferrnal"
            )
        ),
        ref
    )
})

def <- StockholmAAMultipleAlignment(
    Biostrings::AAMultipleAlignment(
        c(
            "O83071/192-246" = "MTCRAQLIAVPRASSLAE..AIACAQKM....RVSRVPVYERS",
            "O83071/259-312" = "MQHVSAPVFVFECTRLAY..VQHKLRAH....SRAVAIVLDEY",
            "O31698/18-71"   = "MIEADKVAHVQVGNNLEH..ALLVLTKT....GYTAIPVLDPS",
            "O31698/88-139"  = "EVMLTDIPRLHINDPIMK..GFGMVINN......GFVCVENDE",
            "O31699/88-139"  = "EVMLTDIPRLHINDPIMK..GFGMVINN......GFVCVENDE"
        )
    ),
    GF = c(
        "ID" = "CBS",
        "AC" = "PF00571",
        "DE" = "CBS domain",
        "AU" = "Bateman A",
        "CC" = paste("CBS domains are small intracellular modules mostly found",
                     "in 2 or four copies within a protein. "),
        "SQ" = "67"
    ),
    GS = list(
        AC = c(
            "O31698/18-71" = "O31698",
            "O83071/192-246" = "O83071",
            "O83071/259-312" = "O83071",
            "O31698/88-139" = "O31698"
        ),
        OS = c(
            "O31698/88-139" = "Bacillus subtilis"
        )
    ),
    GR = list(
        SA = c(
            "O83071/192-246" = "999887756453524252..55152525....36463774777"
        ),
        SS = c(
            "O83071/259-312" = "CCCCCHHHHHHHHHHHHH..EEEEEEEE....EEEEEEEEEEE",
            "O31698/18-71"   = "CCCHHHHHHHHHHHHHHH..EEEEEEEE....EEEEEEEEHHH",
            "O31698/88-139"  = "CCCCCCCHHHHHHHHHHH..HEEEEEEE....EEEEEEEEEEH"
        ),
        AS = c(
            "O31699/88-139"  = "________________*__________________________"
        ),
        IN = c(
            "O31699/88-139"  = "____________1______________2__________0____"
        )
    ),
    GC = c(
        "SS_cons"            = "CCCCCHHHHHHHHHHHHH..EEEEEEEE....EEEEEEEEEEH"
    )
)

test_that("reads example from format definition", {
  expect_stockholm_equal(read_stockholm_msa("test.stk", type = "AA"), def)
})
