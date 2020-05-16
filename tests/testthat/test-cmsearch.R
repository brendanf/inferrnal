test_that("cmsearch works", {
    expect_known_value(
        cmsearch(cm = cm_5_8S(), seq = sample_rRNA_fasta(), cpu = 1, quiet = TRUE),
        file = "cmsearch"
    )
})

test_that("cmsearch is quiet in quiet mode", {
    expect_silent(
        cmsearch(cm = cm_5_8S(), seq = sample_rRNA_fasta(), cpu = 1, quiet = TRUE)
    )
})

dnass <- Biostrings::readDNAStringSet(sample_rRNA_fasta())
rnass <- Biostrings::RNAStringSet(dnass)
rnachar <- as.character(rnass)
dnachar <- as.character(dnass)
dnasr <- ShortRead::ShortRead(dnass, Biostrings::BStringSet(names(dnass)))
test_that("cmsearch works for different input formats", {
    expect_known_value(cmsearch(cm_5_8S(), rnass, cpu = 1),
                       file = "cmsearch")
    expect_known_value(cmsearch(cm_5_8S(), dnass, cpu = 1),
                       file = "cmsearch")
    expect_known_value(cmsearch(cm_5_8S(), rnachar, cpu = 1),
                       file = "cmsearch")
    expect_known_value(cmsearch(cm_5_8S(), dnachar, cpu = 1),
                       file = "cmsearch")
    expect_known_value(cmsearch(cm_5_8S(), dnasr, cpu = 1),
                       file = "cmsearch")
    expect_error(cmsearch(cm_5_8S(), "bad alphabet", cpu = 1),
                       file = "cmsearch")
    expect_error(cmsearch(cm_5_8S(), 15.3, cpu = 1),
                 file = "cmsearch")
})
