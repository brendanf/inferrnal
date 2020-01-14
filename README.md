
<!-- README.md is generated from README.Rmd. Please edit that file -->

# inferrnal

<!-- badges: start -->

<!-- badges: end -->

R Interface to Call Programs from Infernal RNA Covariance Model Package

Covariance Models (CM) are stochastic models of RNA sequence and
secondary structure. [Infernal](http://eddylab.org/infernal/) (INFERence
of RNA ALignment) is a software package with various command-line tools
related to CMs. `inferrnal` (with two “`r`”s) is a lightweight R
interface which calls the Infernal tools and imports the results to R.
It is developed independently from Infernal, and Infernal must be
installed in order for it to function.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("brendanf/inferrnal")
```

## Examples

So far only two of the tools are implemented.

### cmsearch

In order to search, we need a CM. [Rfam](https://rfam.xfam.org/) has a
wide variety. For this example, we will use the eukaryotic 5.8S rRNA,
with the Rfam ID [`RF00002`](https://rfam.xfam.org/family/RF00002). The
CM is [available from Rfam](https://rfam.xfam.org/family/RF00002/cm),
but it is also included as example data in `inferrnal`.

``` r
cm5_8S <- system.file(file.path("extdata", "RF00002.cm"), package = "inferrnal")
```

We also need some sequences to search. The sample data is from a soil
metabarcoding study focused on fungi. The targetted region includes 5.8S
as well as some of the surrounding rDNA regions.

``` r
sampfasta <- system.file(file.path("extdata", "sample.fasta"), package = "inferrnal")
```

Use `cmsearch()` to locate the 5.8S RNA in each sequence.

``` r
library(inferrnal)
cmsearch(cm = cm5_8S, seq = sampfasta, cpu = 1)
#> # A tibble: 438 x 18
#>    target_name target_accession query_name query_accession mdl   mdl_from mdl_to
#>    <chr>       <chr>            <chr>      <chr>           <chr>    <int>  <int>
#>  1 seq61       -                5_8S_rRNA  RF00002         cm           1    154
#>  2 seq45       -                5_8S_rRNA  RF00002         cm           1    154
#>  3 seq175      -                5_8S_rRNA  RF00002         cm           1    154
#>  4 seq3        -                5_8S_rRNA  RF00002         cm           1    154
#>  5 seq92       -                5_8S_rRNA  RF00002         cm           1    154
#>  6 seq2        -                5_8S_rRNA  RF00002         cm           1    154
#>  7 seq28       -                5_8S_rRNA  RF00002         cm           1    154
#>  8 seq63       -                5_8S_rRNA  RF00002         cm           1    154
#>  9 seq83       -                5_8S_rRNA  RF00002         cm           1    154
#> 10 seq154      -                5_8S_rRNA  RF00002         cm           1    154
#> # … with 428 more rows, and 11 more variables: seq_from <int>, seq_to <int>,
#> #   strand <chr>, trunc <chr>, pass <int>, gc <dbl>, bias <dbl>, score <dbl>,
#> #   E_value <dbl>, inc <chr>, description <chr>
```

Instead of passing a file name, you can also supply a `DNAStringSet` or
`RNAStringSet` object from the `Biostrings` package.

``` r
sampseqs <- Biostrings::readDNAStringSet(sampfasta)
cmsearch(cm = cm5_8S, seq = sampseqs, cpu = 1)
#> # A tibble: 438 x 18
#>    target_name target_accession query_name query_accession mdl   mdl_from mdl_to
#>    <chr>       <chr>            <chr>      <chr>           <chr>    <int>  <int>
#>  1 seq61       -                5_8S_rRNA  RF00002         cm           1    154
#>  2 seq45       -                5_8S_rRNA  RF00002         cm           1    154
#>  3 seq175      -                5_8S_rRNA  RF00002         cm           1    154
#>  4 seq3        -                5_8S_rRNA  RF00002         cm           1    154
#>  5 seq92       -                5_8S_rRNA  RF00002         cm           1    154
#>  6 seq2        -                5_8S_rRNA  RF00002         cm           1    154
#>  7 seq28       -                5_8S_rRNA  RF00002         cm           1    154
#>  8 seq63       -                5_8S_rRNA  RF00002         cm           1    154
#>  9 seq83       -                5_8S_rRNA  RF00002         cm           1    154
#> 10 seq154      -                5_8S_rRNA  RF00002         cm           1    154
#> # … with 428 more rows, and 11 more variables: seq_from <int>, seq_to <int>,
#> #   strand <chr>, trunc <chr>, pass <int>, gc <dbl>, bias <dbl>, score <dbl>,
#> #   E_value <dbl>, inc <chr>, description <chr>
```

`cmsearch`, by default, returns a table with information about each hit.
However, it can optionally also output an alignment of the hits in
Stockholm format.

``` r
alnfile <- tempfile("alignment-", fileext = ".stk")
cmsearch(cm = cm5_8S, seq = sampseqs, alignment = alnfile)
#> # A tibble: 438 x 18
#>    target_name target_accession query_name query_accession mdl   mdl_from mdl_to
#>    <chr>       <chr>            <chr>      <chr>           <chr>    <int>  <int>
#>  1 seq61       -                5_8S_rRNA  RF00002         cm           1    154
#>  2 seq45       -                5_8S_rRNA  RF00002         cm           1    154
#>  3 seq175      -                5_8S_rRNA  RF00002         cm           1    154
#>  4 seq3        -                5_8S_rRNA  RF00002         cm           1    154
#>  5 seq92       -                5_8S_rRNA  RF00002         cm           1    154
#>  6 seq2        -                5_8S_rRNA  RF00002         cm           1    154
#>  7 seq28       -                5_8S_rRNA  RF00002         cm           1    154
#>  8 seq63       -                5_8S_rRNA  RF00002         cm           1    154
#>  9 seq83       -                5_8S_rRNA  RF00002         cm           1    154
#> 10 seq154      -                5_8S_rRNA  RF00002         cm           1    154
#> # … with 428 more rows, and 11 more variables: seq_from <int>, seq_to <int>,
#> #   strand <chr>, trunc <chr>, pass <int>, gc <dbl>, bias <dbl>, score <dbl>,
#> #   E_value <dbl>, inc <chr>, description <chr>
```

`inferrnal` includes a simple parser for Stockholm alignments, which
also imports column annotations.

``` r
msa <- read_stockholm_msa(alnfile)
```

`read_stockholm_msa` returns a named list. The alignment itself is given
as an `RNAMultipleSequenceAlignment` object in element `alignment`.

``` r
msa$alignment
#> RNAMultipleAlignment with 438 rows and 235 columns
#>        aln                                                  names               
#>   [1] AACUUUC.A.GCAA.CGG.AU..C.....UGCC...UGU.UUGA.G.UGUCG seq61/192-345
#>   [2] AACUUUC.A.GCAA.CGG.AU..C.....UGCC...UGC.UUGA.G.UGUCA seq45/295-448
#>   [3] AACUUUC.A.GCAA.CGG.AU..C.....UGCC...UGU.UUGA.G.UGUCA seq175/292-445
#>   [4] AACUUUC.A.ACAA.CGG.AU..C.....UGCC...UGU.UUGA.G.UGUCA seq3/236-389
#>   [5] AACUUUC.A.GCAA.CGG.AU..C.....UGCC...UGU.UUGA.G.UGUCA seq92/209-362
#>   [6] AACUUUC.A.GCAA.CGG.AU..C.....UGCC...UGU.UUGA.G.UGUCG seq2/193-346
#>   [7] AACUUUC.A.GCAA.CGG.AU..C.....UGCC...UGU.UUGA.G.UGUCG seq28/192-345
#>   [8] AACUUUC.A.GCAA.CGG.AU..C.....UGCC...UGU.UUGA.G.UGUCG seq63/192-345
#>   [9] AACUUUC.A.GCAA.CGG.AU..C.....UGCC...UGU.UUGA.G.UGUCG seq83/192-345
#>   ... ...
#> [430] AACACGC.A.ACGG.UGG.AC..C.....UUGC...UUG.UUGA.G.CCUGG seq43/218-373
#> [431] AACUGUU.G.GCAA.UGG.AU..C.....CAUC...AUU.GUCA.G.UGCUG seq76/126-272
#> [432] AACCCUU.A.ACGG.UGG.AU..G.....UACU...UGU.UUGA.G.CUCGG seq424/124-285
#> [433] AACUGUU.G.GCAA.CGG.AU..C.....CAUC...AUU.GUCA.G.UGCUG seq211/127-273
#> [434] AACGCAU.A.UCCG.UGG.AU..G.....UACC...UGG.CUAA.G.CUCGG seq364/100-257
#> [435] AACUUCA.A.CACA.CGG.AU..U.....CAAC..CUGU.UCGA.G.CGCCA seq452/179-332
#> [436] AAGGCAU.A.ACGG.UGG.AU..G.....UGUU...CUG.CUGA.G.CUGUA seq174/171-324
#> [437] AACUCUG.U.GCAA.UGG.AG..C.....CUUG...CGU.UAGA.G.UGCUA seq325/57-201
#> [438] UUGGCUC.A.GCGG.UGG.AU..A.....UUCU...UGU.GUCA.A.AGCUU seq276/85-237
```

Other useful information which Infernal outputs are the consensus
secondary structure and the reference annotation, both of which are
defined by the CM:

``` r
msa$SS_cons
#> [1] ":::::::.:.::::.:::.::..:.:.:::.:::...:::.:.::.::.::::.::.((((.<<<.<.__.____.>>.>>,.,.,.,.,,.,.<<<-..---.<<._.__._..>>.--.-----.-.-.-.>>..>.,,,,,.,)))).,,.,<<.<.___>>>.<<<<.<..<<<<.........____......>>>>>.>.>>.>:.::::...:::.::::.:.:::::"
msa$RF
#> [1] "AACuuUu.A.gCGA.UGG.AU..g.u.CUu.GGC...UCc.c.Gu.aU.CGAU.GA.Agaa.CGC.a.GC.aAAa.uG.CGA.U.A.c.GU.a.guGU..GAA.uu.G.CA.G..aa.Uu.ccgUg.A.A.U.Ca..c.CGAAu.cuucG.AA.CGC.a.aaUuGC.Gccc.c..cggg.........Uuuu......cccgg.g.gg.CA.Ugcc...UGu.uugA.G.UGUCa"
```

### cmalign

If you have sequences which have already been trimmed to contain only
the RNA defined by the CM (possibly truncated, but not extended), then
you can align them to the CM using `cmalign`. This is much faster than
`cmsearch`. This example uses the results of `cmsearch` from the
previous section, after removing gaps.

``` r
unaln <- system.file(file.path("extdata", "samp_5_8S.fasta"), package = "inferrnal")
unaln_seq <- Biostrings::readRNAStringSet(unaln)
unaln_seq
#>   A RNAStringSet instance of length 438
#>       width seq                                             names               
#>   [1]   154 AACUUUCAGCAACGGAUCUCUU...GGGCAUGCCUGUUUGAGUGUCG seq61/192-345
#>   [2]   154 AACUUUCAGCAACGGAUCUCUU...GAGCAUGCCUGCUUGAGUGUCA seq45/295-448
#>   [3]   154 AACUUUCAGCAACGGAUCUCUU...GAGCAUGCCUGUUUGAGUGUCA seq175/292-445
#>   [4]   154 AACUUUCAACAACGGAUCUCUU...GAGCAUGCCUGUUUGAGUGUCA seq3/236-389
#>   [5]   154 AACUUUCAGCAACGGAUCUCUU...GAGCAUGCCUGUUUGAGUGUCA seq92/209-362
#>   ...   ... ...
#> [434]   158 AACGCAUAUCCGUGGAUGACUC...GGCCGUACCUGGCUAAGCUCGG seq364/100-257
#> [435]   154 AACUUCAACACACGGAUUCUUG...GGCACAACCUGUUCGAGCGCCA seq452/179-332
#> [436]   154 AAGGCAUAACGGUGGAUGAAUU...AGGUGUGUUCUGCUGAGCUGUA seq174/171-324
#> [437]   145 AACUCUGUGCAAUGGAGCACAC...AGGUACUUGCGUUAGAGUGCUA seq325/57-201
#> [438]   153 UUGGCUCAGCGGUGGAUAUCUA...GUGAAUUCUUGUGUCAAAGCUU seq276/85-237
aln <- cmalign(cm5_8S, unaln_seq, cpu = 1)
aln$alignment
#> RNAMultipleAlignment with 438 rows and 235 columns
#>        aln                                                  names               
#>   [1] AACUUUC.A.GCAA.CGG.AU..C.....UGCC...UGU.UUGA.G.UGUCG seq61/192-345
#>   [2] AACUUUC.A.GCAA.CGG.AU..C.....UGCC...UGC.UUGA.G.UGUCA seq45/295-448
#>   [3] AACUUUC.A.GCAA.CGG.AU..C.....UGCC...UGU.UUGA.G.UGUCA seq175/292-445
#>   [4] AACUUUC.A.ACAA.CGG.AU..C.....UGCC...UGU.UUGA.G.UGUCA seq3/236-389
#>   [5] AACUUUC.A.GCAA.CGG.AU..C.....UGCC...UGU.UUGA.G.UGUCA seq92/209-362
#>   [6] AACUUUC.A.GCAA.CGG.AU..C.....UGCC...UGU.UUGA.G.UGUCG seq2/193-346
#>   [7] AACUUUC.A.GCAA.CGG.AU..C.....UGCC...UGU.UUGA.G.UGUCG seq28/192-345
#>   [8] AACUUUC.A.GCAA.CGG.AU..C.....UGCC...UGU.UUGA.G.UGUCG seq63/192-345
#>   [9] AACUUUC.A.GCAA.CGG.AU..C.....UGCC...UGU.UUGA.G.UGUCG seq83/192-345
#>   ... ...
#> [430] AACACGC.A.ACGG.UGG.AC..C.....UUGC...UUG.UUGA.G.CCUGG seq43/218-373
#> [431] AACUGUU.G.GCAA.UGG.AU..C.....CAUC...AUU.GUCA.G.UGCUG seq76/126-272
#> [432] AACCCUU.A.ACGG.UGG.AU..G.....UACU...UGU.UUGA.G.CUCGG seq424/124-285
#> [433] AACUGUU.G.GCAA.CGG.AU..C.....CAUC...AUU.GUCA.G.UGCUG seq211/127-273
#> [434] AACGCAU.A.UCCG.UGG.AU..G.....UACC...UGG.CUAA.G.CUCGG seq364/100-257
#> [435] AACUUCA.A.CACA.CGG.AU..U.....CAAC..CUGU.UCGA.G.CGCCA seq452/179-332
#> [436] AAGGCAU.A.ACGG.UGG.AU..G.....UGUU...CUG.CUGA.G.CUGUA seq174/171-324
#> [437] AACUCUG.U.GCAA.UGG.AG..C.....CUUG...CGU.UAGA.G.UGCUA seq325/57-201
#> [438] UUGGCUC.A.GCGG.UGG.AU..A.....UUCU...UGU.GUCA.A.AGCUU seq276/85-237
aln$SS_cons
#> [1] ":::::::.:.::::.:::.::..:.:.:::.:::...:::.:.::.::.::::.::.((((.<<<.<.__.____.>>.>>,.,.,.,.,,.,.<<<-..---.<<._.__._..>>.--.-----.-.-.-.>>..>.,,,,,.,)))).,,.,<<.<.___>>>.<<<<.<..<<<<.........____......>>>>>.>.>>.>:.::::...:::.::::.:.:::::"
aln$RF
#> [1] "AACuuUu.A.gCGA.UGG.AU..g.u.CUu.GGC...UCc.c.Gu.aU.CGAU.GA.Agaa.CGC.a.GC.aAAa.uG.CGA.U.A.c.GU.a.guGU..GAA.uu.G.CA.G..aa.Uu.ccgUg.A.A.U.Ca..c.CGAAu.cuucG.AA.CGC.a.aaUuGC.Gccc.c..cggg.........Uuuu......cccgg.g.gg.CA.Ugcc...UGu.uugA.G.UGUCa"
```
