
<!-- README.md is generated from README.Rmd. Please edit that file -->

# inferrnal

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/brendanf/inferrnal.svg?branch=master)](https://travis-ci.com/brendanf/inferrnal)
[![Codecov test
coverage](https://codecov.io/gh/brendanf/inferrnal/branch/master/graph/badge.svg)](https://codecov.io/gh/brendanf/inferrnal?branch=master)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- badges: end -->

R Interface to Call Programs from Infernal RNA Covariance Model Package

Covariance Models (CM) are stochastic models of RNA sequence and
secondary structure. [Infernal](http://eddylab.org/infernal/) (INFERence
of RNA ALignment) is a software package with various command-line tools
related to CMs. `inferrnal` (with two “`r`”s) is a lightweight R
interface which calls the Infernal tools and imports the results to R.
It is developed independently from Infernal, and Infernal must be
installed in order for it to function. Note that Infernal does not work
on Windows.

## Installation

After installing Infernal, you can install the development version of
`inferrnal` from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("brendanf/inferrnal")
```

Submission to Bioconductor is pending.

## Examples

So far only two of the tools are implemented.

### cmsearch

In order to search, we need a CM. [Rfam](https://rfam.xfam.org/) has a
wide variety. For this example, we will use the eukaryotic 5.8S rRNA,
with the Rfam ID [`RF00002`](https://rfam.xfam.org/family/RF00002). The
CM is [available from Rfam](https://rfam.xfam.org/family/RF00002/cm),
but it is also included as example data in `inferrnal`.

``` r
library(inferrnal)
cm <- cm_5_8S()
```

We also need some sequences to search. The sample data is from a soil
metabarcoding study focused on fungi. The targeted region includes 5.8S
as well as some of the surrounding rDNA regions.

``` r
sampfasta <- sample_rRNA_fasta()
```

Use `cmsearch()` to locate the 5.8S RNA in each sequence.

``` r
library(inferrnal)
cmsearch(cm = cm, seq = sampfasta, cpu = 1)
#> # A tibble: 49 x 18
#>    target_name target_accession query_name query_accession mdl   mdl_from mdl_to
#>    <chr>       <chr>            <chr>      <chr>           <chr>    <int>  <int>
#>  1 seq45       -                5_8S_rRNA  RF00002         cm           1    154
#>  2 seq3        -                5_8S_rRNA  RF00002         cm           1    154
#>  3 seq2        -                5_8S_rRNA  RF00002         cm           1    154
#>  4 seq28       -                5_8S_rRNA  RF00002         cm           1    154
#>  5 seq23       -                5_8S_rRNA  RF00002         cm           1    154
#>  6 seq9        -                5_8S_rRNA  RF00002         cm           1    154
#>  7 seq7        -                5_8S_rRNA  RF00002         cm           1    154
#>  8 seq48       -                5_8S_rRNA  RF00002         cm           1    154
#>  9 seq5        -                5_8S_rRNA  RF00002         cm           1    154
#> 10 seq6        -                5_8S_rRNA  RF00002         cm           1    154
#> # … with 39 more rows, and 11 more variables: seq_from <int>, seq_to <int>,
#> #   strand <chr>, trunc <chr>, pass <int>, gc <dbl>, bias <dbl>, score <dbl>,
#> #   E_value <dbl>, inc <chr>, description <chr>
```

Instead of passing a file name, you can also supply a `DNAStringSet` or
`RNAStringSet` object from the `Biostrings` package.

``` r
sampseqs <- Biostrings::readDNAStringSet(sampfasta)
cmsearch(cm = cm, seq = sampseqs, cpu = 1)
#> # A tibble: 49 x 18
#>    target_name target_accession query_name query_accession mdl   mdl_from mdl_to
#>    <chr>       <chr>            <chr>      <chr>           <chr>    <int>  <int>
#>  1 seq45       -                5_8S_rRNA  RF00002         cm           1    154
#>  2 seq3        -                5_8S_rRNA  RF00002         cm           1    154
#>  3 seq2        -                5_8S_rRNA  RF00002         cm           1    154
#>  4 seq28       -                5_8S_rRNA  RF00002         cm           1    154
#>  5 seq23       -                5_8S_rRNA  RF00002         cm           1    154
#>  6 seq9        -                5_8S_rRNA  RF00002         cm           1    154
#>  7 seq7        -                5_8S_rRNA  RF00002         cm           1    154
#>  8 seq48       -                5_8S_rRNA  RF00002         cm           1    154
#>  9 seq5        -                5_8S_rRNA  RF00002         cm           1    154
#> 10 seq6        -                5_8S_rRNA  RF00002         cm           1    154
#> # … with 39 more rows, and 11 more variables: seq_from <int>, seq_to <int>,
#> #   strand <chr>, trunc <chr>, pass <int>, gc <dbl>, bias <dbl>, score <dbl>,
#> #   E_value <dbl>, inc <chr>, description <chr>
```

`cmsearch`, by default, returns a table with information about each hit.
However, it can optionally also output an alignment of the hits in
Stockholm format.

``` r
alnfile <- tempfile("alignment-", fileext = ".stk")
cmsearch(cm = cm, seq = sampseqs, alignment = alnfile)
#> # A tibble: 49 x 18
#>    target_name target_accession query_name query_accession mdl   mdl_from mdl_to
#>    <chr>       <chr>            <chr>      <chr>           <chr>    <int>  <int>
#>  1 seq45       -                5_8S_rRNA  RF00002         cm           1    154
#>  2 seq3        -                5_8S_rRNA  RF00002         cm           1    154
#>  3 seq2        -                5_8S_rRNA  RF00002         cm           1    154
#>  4 seq28       -                5_8S_rRNA  RF00002         cm           1    154
#>  5 seq23       -                5_8S_rRNA  RF00002         cm           1    154
#>  6 seq9        -                5_8S_rRNA  RF00002         cm           1    154
#>  7 seq7        -                5_8S_rRNA  RF00002         cm           1    154
#>  8 seq48       -                5_8S_rRNA  RF00002         cm           1    154
#>  9 seq5        -                5_8S_rRNA  RF00002         cm           1    154
#> 10 seq6        -                5_8S_rRNA  RF00002         cm           1    154
#> # … with 39 more rows, and 11 more variables: seq_from <int>, seq_to <int>,
#> #   strand <chr>, trunc <chr>, pass <int>, gc <dbl>, bias <dbl>, score <dbl>,
#> #   E_value <dbl>, inc <chr>, description <chr>
```

`inferrnal` includes a simple parser for Stockholm alignments, which
also imports column annotations.

``` r
msa <- read_stockholm_msa(alnfile)
```

`read_stockholm_msa` returns a named list. The alignment itself is given
as an `RNAMultipleAlignment` object in element `alignment`.

``` r
msa$alignment
#> RNAMultipleAlignment with 48 rows and 164 columns
#>       aln                                                   names               
#>  [1] AACUUUCAGCAACGGAUCUCUUGGC...GAGGAGCAUGCCUGCUUGAGUGUCA seq45/295-448
#>  [2] AACUUUCAACAACGGAUCUCUUGGC...GAGGAGCAUGCCUGUUUGAGUGUCA seq3/236-389
#>  [3] AACUUUCAGCAACGGAUCUCUUGGC...GGAGGGCAUGCCUGUUUGAGUGUCG seq2/193-346
#>  [4] AACUUUCAGCAACGGAUCUCUUGGC...GGAGGGCAUGCCUGUUUGAGUGUCG seq28/192-345
#>  [5] AACUUUCAACAACGGAUCUCUUGGU...GAGGGGCAUGCCUGUUCGAGCGUCA seq23/194-347
#>  [6] AACUUUCAACAACGGAUCUCUUGGU...GAGGGGCAUGCCUGUUCGAGCGUCA seq9/170-323
#>  [7] AACUUUCAACAAUGGAUCUCUUGGU...GGGGGGCAUGCCUGUCCGAGCGUCA seq7/191-344
#>  [8] AACUUUCAGCAACGGAUCUCUUGGC...GGAGGGCAUGCCUGUUUGAGUGUCG seq48/192-345
#>  [9] AACUUUCAGCAACGGAUCUCUUGGC...GAGGAGUAUGCCUGUUUGAGUAUCA seq5/257-410
#>  ... ...
#> [40] AACUUUCAGCAACGGAUCUCUUGGC...GAGGAGUAUGCCUGUUUGAGUAUCA seq41/258-414
#> [41] AACUUUCAGCAACGGAUCUCUUGGC...UGUGAGUACACUUGUUUGAGCGUCA seq29/319-475
#> [42] CACUAUUAGCGAUGGAUGUCUUGGA...CAGCAGUAGGUUGGUCUCAGCAUCU seq49/188-340
#> [43] GACUCCCGGCAACGGAUAUCUCGGC...CGAGGGCACGCCUGCCUCUUGGGCG seq20/234-388
#> [44] AACUCUCAGCGAUGGAUGACUCGAC...CUGAAGUAUGUUUGGCUCGGUAUCA seq4/95-248
#> [45] AACUCUCAGCGAUGGAUGACUCGAC...CUGAAGUAUGUUUGGCUCGGUAUCA seq32/94-246
#> [46] AACUCUCAGCGAUGGAUGACUCGAC...CUGAAGUAUGUUUGGCUCGGUAUCA seq24/95-246
#> [47] UAGCAUCAGCGAUUAACGUCUUGGU...AUUGAGUGCACUUGCUUCAGUGUGG seq37/93-246
#> [48] AACACGCAACGGUGGACCACUCGGC...GCCAGCUCUUGCUUGUUGAGCCUGG seq43/218-373
```

Other useful information which Infernal outputs are the consensus
secondary structure and the reference annotation, both of which are
defined by the CM:

``` r
msa$SS_cons
#> [1] ":::::::::::::::::::::::::::::::::::::::((((<<<<__.____.>>.>>,,,,.,,,.<<<-..---<<____>>---------->>>,,,,,,)))),,,<<<___>>><<<<<<<<<..____.>>>>>>>>>::::::::::::::::::"
msa$RF
#> [1] "AACuuUuAgCGAUGGAUguCUuGGCUCccGuaUCGAUGAAgaaCGCaGC.aAAa.uG.CGAUAc.GUa.guGU..GAAuuGCAGaaUuccgUgAAUCacCGAAucuucGAACGCaaaUuGCGcccccggg..Uuuu.cccgggggCAUgccUGuuugAGUGUCa"
```

### cmalign

If you have sequences which have already been trimmed to contain only
the RNA defined by the CM (possibly truncated, but not extended), then
you can align them to the CM using `cmalign`. This is much faster than
`cmsearch`. This example uses the results of `cmsearch` from the
previous section, after removing gaps.

``` r
unaln <- sample_rRNA_5_8S()
unaln_seq <- Biostrings::readRNAStringSet(unaln)
unaln_seq
#>   A RNAStringSet instance of length 48
#>      width seq                                              names               
#>  [1]   154 AACUUUCAGCAACGGAUCUCUUG...GAGCAUGCCUGCUUGAGUGUCA seq45/295-448
#>  [2]   154 AACUUUCAACAACGGAUCUCUUG...GAGCAUGCCUGUUUGAGUGUCA seq3/236-389
#>  [3]   154 AACUUUCAGCAACGGAUCUCUUG...GGGCAUGCCUGUUUGAGUGUCG seq2/193-346
#>  [4]   154 AACUUUCAGCAACGGAUCUCUUG...GGGCAUGCCUGUUUGAGUGUCG seq28/192-345
#>  [5]   154 AACUUUCAACAACGGAUCUCUUG...GGGCAUGCCUGUUCGAGCGUCA seq23/194-347
#>  ...   ... ...
#> [44]   154 AACUCUCAGCGAUGGAUGACUCG...AAGUAUGUUUGGCUCGGUAUCA seq4/95-248
#> [45]   153 AACUCUCAGCGAUGGAUGACUCG...AAGUAUGUUUGGCUCGGUAUCA seq32/94-246
#> [46]   152 AACUCUCAGCGAUGGAUGACUCG...AAGUAUGUUUGGCUCGGUAUCA seq24/95-246
#> [47]   154 UAGCAUCAGCGAUUAACGUCUUG...GAGUGCACUUGCUUCAGUGUGG seq37/93-246
#> [48]   156 AACACGCAACGGUGGACCACUCG...AGCUCUUGCUUGUUGAGCCUGG seq43/218-373
aln <- cmalign(cm, unaln_seq, cpu = 1)
aln$alignment
#> RNAMultipleAlignment with 48 rows and 164 columns
#>       aln                                                   names               
#>  [1] AACUUUCAGCAACGGAUCUCUUGGC...GAGGAGCAUGCCUGCUUGAGUGUCA seq45/295-448
#>  [2] AACUUUCAACAACGGAUCUCUUGGC...GAGGAGCAUGCCUGUUUGAGUGUCA seq3/236-389
#>  [3] AACUUUCAGCAACGGAUCUCUUGGC...GGAGGGCAUGCCUGUUUGAGUGUCG seq2/193-346
#>  [4] AACUUUCAGCAACGGAUCUCUUGGC...GGAGGGCAUGCCUGUUUGAGUGUCG seq28/192-345
#>  [5] AACUUUCAACAACGGAUCUCUUGGU...GAGGGGCAUGCCUGUUCGAGCGUCA seq23/194-347
#>  [6] AACUUUCAACAACGGAUCUCUUGGU...GAGGGGCAUGCCUGUUCGAGCGUCA seq9/170-323
#>  [7] AACUUUCAACAAUGGAUCUCUUGGU...GGGGGGCAUGCCUGUCCGAGCGUCA seq7/191-344
#>  [8] AACUUUCAGCAACGGAUCUCUUGGC...GGAGGGCAUGCCUGUUUGAGUGUCG seq48/192-345
#>  [9] AACUUUCAGCAACGGAUCUCUUGGC...GAGGAGUAUGCCUGUUUGAGUAUCA seq5/257-410
#>  ... ...
#> [40] AACUUUCAGCAACGGAUCUCUUGGC...GAGGAGUAUGCCUGUUUGAGUAUCA seq41/258-414
#> [41] AACUUUCAGCAACGGAUCUCUUGGC...UGUGAGUACACUUGUUUGAGCGUCA seq29/319-475
#> [42] CACUAUUAGCGAUGGAUGUCUUGGA...CAGCAGUAGGUUGGUCUCAGCAUCU seq49/188-340
#> [43] GACUCCCGGCAACGGAUAUCUCGGC...CGAGGGCACGCCUGCCUCUUGGGCG seq20/234-388
#> [44] AACUCUCAGCGAUGGAUGACUCGAC...CUGAAGUAUGUUUGGCUCGGUAUCA seq4/95-248
#> [45] AACUCUCAGCGAUGGAUGACUCGAC...CUGAAGUAUGUUUGGCUCGGUAUCA seq32/94-246
#> [46] AACUCUCAGCGAUGGAUGACUCGAC...CUGAAGUAUGUUUGGCUCGGUAUCA seq24/95-246
#> [47] UAGCAUCAGCGAUUAACGUCUUGGU...AUUGAGUGCACUUGCUUCAGUGUGG seq37/93-246
#> [48] AACACGCAACGGUGGACCACUCGGC...GCCAGCUCUUGCUUGUUGAGCCUGG seq43/218-373
aln$SS_cons
#> [1] ":::::::::::::::::::::::::::::::::::::::((((<<<<__.____.>>.>>,,,,.,,,.<<<-..---<<____>>---------->>>,,,,,,)))),,,<<<___>>><<<<<<<<<..____.>>>>>>>>>::::::::::::::::::"
aln$RF
#> [1] "AACuuUuAgCGAUGGAUguCUuGGCUCccGuaUCGAUGAAgaaCGCaGC.aAAa.uG.CGAUAc.GUa.guGU..GAAuuGCAGaaUuccgUgAAUCacCGAAucuucGAACGCaaaUuGCGcccccggg..Uuuu.cccgggggCAUgccUGuuugAGUGUCa"
```

## cmbuild

`cmbuild`, included in `inferrnal` since version 0.99.5, is used to
create new CMs from annotated multiple sequence alignments. To
illustrate the process, we use the seed alignment for the 5.8S rRNA CM
from RFAM. It is included as a sample file in `inferrnal`.

``` r
new_cm <- file.path(tempdir(), "5_8S.cm")
cmbuild(new_cm, msafile = stk_5_8S(), force = TRUE, quiet = FALSE)
```

This CM is not calibrated, so it cannot be used for `cmsearch`, but it
can be used in `cmalign`.

``` r
aln2 <- cmalign(new_cm, unaln_seq, cpu = 1)
```

The resulting alignment is the same as the one created using the CM from
RFAM, because they are based on the same seed alignment, and used the
same (default) options for `cmbuild`.

``` r
all.equal(aln, aln2)
#> [1] TRUE
```
