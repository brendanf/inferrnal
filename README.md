
<!-- README.md is generated from README.Rmd. Please edit that file -->

# inferrnal

<!-- badges: start -->

[![R build
status](https://github.com/brendanf/inferrnal/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/brendanf/inferrnal/actions)
[![Codecov test
coverage](https://codecov.io/gh/brendanf/inferrnal/branch/master/graph/badge.svg)](https://codecov.io/gh/brendanf/inferrnal?branch=master)
[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
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

### Installing Infernal

The required Infernal package can be installed from Bioconda:

    conda install -c bioconda infernal

or in Debian/Ubuntu Linux using apt:

    sudo apt-get install infernal

or using Homebrew in MacOS:

    brew tap brewsci/bio
    brew install infernal

For other installation options, including source installation, see the
[Infernal homepage](http://eddylab.org/infernal/).

### Installing `inferrnal`

After Infernal is installed, the released version of `inferrnal` can be
installed from Bioconductor:

``` r
if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("inferrnal")
```

Alternatively, the development version of `inferrnal` can be installed
from [GitHub](https://github.com/) with:

``` r
if(!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("brendanf/inferrnal")
```

## Examples

So far three of the tools are implemented: `cmsearch`, `cmalign`, and
`cmbuild`.

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
#>    target_name target_accession query_name query_accession mdl mdl_from mdl_to
#> 1        seq45                -  5_8S_rRNA         RF00002  cm        1    154
#> 2         seq3                -  5_8S_rRNA         RF00002  cm        1    154
#> 3         seq2                -  5_8S_rRNA         RF00002  cm        1    154
#> 4        seq28                -  5_8S_rRNA         RF00002  cm        1    154
#> 5        seq23                -  5_8S_rRNA         RF00002  cm        1    154
#> 6         seq9                -  5_8S_rRNA         RF00002  cm        1    154
#> 7         seq7                -  5_8S_rRNA         RF00002  cm        1    154
#> 8        seq48                -  5_8S_rRNA         RF00002  cm        1    154
#> 9         seq5                -  5_8S_rRNA         RF00002  cm        1    154
#> 10        seq6                -  5_8S_rRNA         RF00002  cm        1    154
#> 11       seq11                -  5_8S_rRNA         RF00002  cm        1    154
#> 12       seq12                -  5_8S_rRNA         RF00002  cm        1    154
#> 13       seq13                -  5_8S_rRNA         RF00002  cm        1    154
#> 14       seq17                -  5_8S_rRNA         RF00002  cm        1    154
#> 15       seq18                -  5_8S_rRNA         RF00002  cm        1    154
#> 16       seq19                -  5_8S_rRNA         RF00002  cm        1    154
#> 17       seq25                -  5_8S_rRNA         RF00002  cm        1    154
#> 18       seq47                -  5_8S_rRNA         RF00002  cm        1    154
#> 19       seq14                -  5_8S_rRNA         RF00002  cm        1    154
#> 20        seq8                -  5_8S_rRNA         RF00002  cm        1    154
#> 21       seq21                -  5_8S_rRNA         RF00002  cm        1    154
#> 22       seq36                -  5_8S_rRNA         RF00002  cm        1    154
#> 23       seq10                -  5_8S_rRNA         RF00002  cm        1    154
#> 24       seq22                -  5_8S_rRNA         RF00002  cm        1    154
#> 25       seq50                -  5_8S_rRNA         RF00002  cm        1    154
#> 26       seq44                -  5_8S_rRNA         RF00002  cm        1    154
#> 27       seq27                -  5_8S_rRNA         RF00002  cm        1    154
#> 28       seq35                -  5_8S_rRNA         RF00002  cm        1    154
#> 29       seq15                -  5_8S_rRNA         RF00002  cm        1    154
#> 30       seq38                -  5_8S_rRNA         RF00002  cm        1    154
#> 31       seq30                -  5_8S_rRNA         RF00002  cm        1    154
#> 32       seq26                -  5_8S_rRNA         RF00002  cm        1    154
#> 33       seq33                -  5_8S_rRNA         RF00002  cm        1    154
#> 34       seq46                -  5_8S_rRNA         RF00002  cm        1    154
#> 35       seq34                -  5_8S_rRNA         RF00002  cm        1    154
#> 36       seq39                -  5_8S_rRNA         RF00002  cm        1    154
#> 37       seq16                -  5_8S_rRNA         RF00002  cm        1    154
#> 38       seq31                -  5_8S_rRNA         RF00002  cm        1    154
#> 39        seq1                -  5_8S_rRNA         RF00002  cm        1    154
#> 40       seq41                -  5_8S_rRNA         RF00002  cm        1    154
#> 41       seq29                -  5_8S_rRNA         RF00002  cm        1    154
#> 42       seq49                -  5_8S_rRNA         RF00002  cm        1    154
#> 43       seq20                -  5_8S_rRNA         RF00002  cm        1    154
#> 44        seq4                -  5_8S_rRNA         RF00002  cm        1    154
#> 45       seq32                -  5_8S_rRNA         RF00002  cm        1    154
#> 46       seq24                -  5_8S_rRNA         RF00002  cm        1    154
#> 47       seq37                -  5_8S_rRNA         RF00002  cm        1    154
#> 48       seq43                -  5_8S_rRNA         RF00002  cm        1    154
#> 49       seq37                -  5_8S_rRNA         RF00002  cm      145    154
#>    seq_from seq_to strand trunc pass   gc bias score E_value inc description
#> 1       295    448      +    no    1 0.47    0 171.3 1.9e-27   !           -
#> 2       236    389      +    no    1 0.45    0 171.0 2.1e-27   !           -
#> 3       193    346      +    no    1 0.48    0 170.7 2.4e-27   !           -
#> 4       192    345      +    no    1 0.48    0 170.7 2.4e-27   !           -
#> 5       194    347      +    no    1 0.47    0 168.5 5.3e-27   !           -
#> 6       170    323      +    no    1 0.47    0 167.3 7.9e-27   !           -
#> 7       191    344      +    no    1 0.50    0 167.3 8.1e-27   !           -
#> 8       192    345      +    no    1 0.49    0 167.1 8.7e-27   !           -
#> 9       257    410      +    no    1 0.43    0 165.5 1.5e-26   !           -
#> 10      257    410      +    no    1 0.43    0 165.5 1.5e-26   !           -
#> 11      258    411      +    no    1 0.43    0 165.5 1.5e-26   !           -
#> 12      257    410      +    no    1 0.43    0 165.5 1.5e-26   !           -
#> 13      256    409      +    no    1 0.43    0 165.5 1.5e-26   !           -
#> 14      257    410      +    no    1 0.43    0 165.5 1.5e-26   !           -
#> 15      258    411      +    no    1 0.43    0 165.5 1.5e-26   !           -
#> 16      257    410      +    no    1 0.43    0 165.5 1.5e-26   !           -
#> 17      256    409      +    no    1 0.43    0 165.5 1.5e-26   !           -
#> 18      257    410      +    no    1 0.43    0 165.5 1.5e-26   !           -
#> 19      157    310      +    no    1 0.46    0 165.3 1.6e-26   !           -
#> 20      154    307      +    no    1 0.47    0 164.7 2.0e-26   !           -
#> 21      256    409      +    no    1 0.42    0 164.7 2.1e-26   !           -
#> 22      256    409      +    no    1 0.44    0 164.1 2.5e-26   !           -
#> 23      257    410      +    no    1 0.44    0 163.9 2.7e-26   !           -
#> 24      256    409      +    no    1 0.44    0 163.9 2.7e-26   !           -
#> 25      257    410      +    no    1 0.44    0 163.7 2.8e-26   !           -
#> 26      255    408      +    no    1 0.44    0 163.6 3.0e-26   !           -
#> 27      216    368      +    no    1 0.46    0 162.5 4.5e-26   !           -
#> 28      178    331      +    no    1 0.47    0 162.4 4.6e-26   !           -
#> 29      258    411      +    no    1 0.42    0 161.0 7.7e-26   !           -
#> 30      156    309      +    no    1 0.47    0 160.6 8.7e-26   !           -
#> 31      195    348      +    no    1 0.48    0 160.1 1.0e-25   !           -
#> 32      188    340      +    no    1 0.50    0 158.7 1.7e-25   !           -
#> 33      223    376      +    no    1 0.50    0 156.6 3.7e-25   !           -
#> 34      230    383      +    no    1 0.49    0 156.6 3.7e-25   !           -
#> 35      219    371      +    no    1 0.45    0 156.5 3.8e-25   !           -
#> 36      263    416      +    no    1 0.44    0 155.7 5.0e-25   !           -
#> 37      256    408      +    no    1 0.43    0 154.0 9.2e-25   !           -
#> 38      194    346      +    no    1 0.48    0 153.9 9.6e-25   !           -
#> 39      115    268      +    no    1 0.45    0 151.2 2.5e-24   !           -
#> 40      258    414      +    no    1 0.44    0 143.4 4.1e-23   !           -
#> 41      319    475      +    no    1 0.45    0 139.6 1.6e-22   !           -
#> 42      188    340      +    no    1 0.46    0 121.1 1.2e-19   !           -
#> 43      234    388      +    no    1 0.57    0 119.6 2.1e-19   !           -
#> 44       95    248      +    no    1 0.51    0 118.1 3.6e-19   !           -
#> 45       94    246      +    no    1 0.50    0 115.1 1.0e-18   !           -
#> 46       95    246      +    no    1 0.50    0 108.2 1.3e-17   !           -
#> 47       93    246      +    no    1 0.44    0 102.4 9.9e-17   !           -
#> 48      218    373      +    no    1 0.51    0  93.5 2.4e-15   !           -
#> 49        1     10      +    5'    2 0.50    0  -5.4 5.9e+00   ?           -
```

Instead of passing a file name, you can also supply a `DNAStringSet` or
`RNAStringSet` object from the `Biostrings` package.

``` r
sampseqs <- Biostrings::readDNAStringSet(sampfasta)
cmsearch(cm = cm, seq = sampseqs, cpu = 1)
#>    target_name target_accession query_name query_accession mdl mdl_from mdl_to
#> 1        seq45                -  5_8S_rRNA         RF00002  cm        1    154
#> 2         seq3                -  5_8S_rRNA         RF00002  cm        1    154
#> 3         seq2                -  5_8S_rRNA         RF00002  cm        1    154
#> 4        seq28                -  5_8S_rRNA         RF00002  cm        1    154
#> 5        seq23                -  5_8S_rRNA         RF00002  cm        1    154
#> 6         seq9                -  5_8S_rRNA         RF00002  cm        1    154
#> 7         seq7                -  5_8S_rRNA         RF00002  cm        1    154
#> 8        seq48                -  5_8S_rRNA         RF00002  cm        1    154
#> 9         seq5                -  5_8S_rRNA         RF00002  cm        1    154
#> 10        seq6                -  5_8S_rRNA         RF00002  cm        1    154
#> 11       seq11                -  5_8S_rRNA         RF00002  cm        1    154
#> 12       seq12                -  5_8S_rRNA         RF00002  cm        1    154
#> 13       seq13                -  5_8S_rRNA         RF00002  cm        1    154
#> 14       seq17                -  5_8S_rRNA         RF00002  cm        1    154
#> 15       seq18                -  5_8S_rRNA         RF00002  cm        1    154
#> 16       seq19                -  5_8S_rRNA         RF00002  cm        1    154
#> 17       seq25                -  5_8S_rRNA         RF00002  cm        1    154
#> 18       seq47                -  5_8S_rRNA         RF00002  cm        1    154
#> 19       seq14                -  5_8S_rRNA         RF00002  cm        1    154
#> 20        seq8                -  5_8S_rRNA         RF00002  cm        1    154
#> 21       seq21                -  5_8S_rRNA         RF00002  cm        1    154
#> 22       seq36                -  5_8S_rRNA         RF00002  cm        1    154
#> 23       seq10                -  5_8S_rRNA         RF00002  cm        1    154
#> 24       seq22                -  5_8S_rRNA         RF00002  cm        1    154
#> 25       seq50                -  5_8S_rRNA         RF00002  cm        1    154
#> 26       seq44                -  5_8S_rRNA         RF00002  cm        1    154
#> 27       seq27                -  5_8S_rRNA         RF00002  cm        1    154
#> 28       seq35                -  5_8S_rRNA         RF00002  cm        1    154
#> 29       seq15                -  5_8S_rRNA         RF00002  cm        1    154
#> 30       seq38                -  5_8S_rRNA         RF00002  cm        1    154
#> 31       seq30                -  5_8S_rRNA         RF00002  cm        1    154
#> 32       seq26                -  5_8S_rRNA         RF00002  cm        1    154
#> 33       seq33                -  5_8S_rRNA         RF00002  cm        1    154
#> 34       seq46                -  5_8S_rRNA         RF00002  cm        1    154
#> 35       seq34                -  5_8S_rRNA         RF00002  cm        1    154
#> 36       seq39                -  5_8S_rRNA         RF00002  cm        1    154
#> 37       seq16                -  5_8S_rRNA         RF00002  cm        1    154
#> 38       seq31                -  5_8S_rRNA         RF00002  cm        1    154
#> 39        seq1                -  5_8S_rRNA         RF00002  cm        1    154
#> 40       seq41                -  5_8S_rRNA         RF00002  cm        1    154
#> 41       seq29                -  5_8S_rRNA         RF00002  cm        1    154
#> 42       seq49                -  5_8S_rRNA         RF00002  cm        1    154
#> 43       seq20                -  5_8S_rRNA         RF00002  cm        1    154
#> 44        seq4                -  5_8S_rRNA         RF00002  cm        1    154
#> 45       seq32                -  5_8S_rRNA         RF00002  cm        1    154
#> 46       seq24                -  5_8S_rRNA         RF00002  cm        1    154
#> 47       seq37                -  5_8S_rRNA         RF00002  cm        1    154
#> 48       seq43                -  5_8S_rRNA         RF00002  cm        1    154
#> 49       seq37                -  5_8S_rRNA         RF00002  cm      145    154
#>    seq_from seq_to strand trunc pass   gc bias score E_value inc description
#> 1       295    448      +    no    1 0.47    0 171.3 1.9e-27   !           -
#> 2       236    389      +    no    1 0.45    0 171.0 2.1e-27   !           -
#> 3       193    346      +    no    1 0.48    0 170.7 2.4e-27   !           -
#> 4       192    345      +    no    1 0.48    0 170.7 2.4e-27   !           -
#> 5       194    347      +    no    1 0.47    0 168.5 5.3e-27   !           -
#> 6       170    323      +    no    1 0.47    0 167.3 7.9e-27   !           -
#> 7       191    344      +    no    1 0.50    0 167.3 8.1e-27   !           -
#> 8       192    345      +    no    1 0.49    0 167.1 8.7e-27   !           -
#> 9       257    410      +    no    1 0.43    0 165.5 1.5e-26   !           -
#> 10      257    410      +    no    1 0.43    0 165.5 1.5e-26   !           -
#> 11      258    411      +    no    1 0.43    0 165.5 1.5e-26   !           -
#> 12      257    410      +    no    1 0.43    0 165.5 1.5e-26   !           -
#> 13      256    409      +    no    1 0.43    0 165.5 1.5e-26   !           -
#> 14      257    410      +    no    1 0.43    0 165.5 1.5e-26   !           -
#> 15      258    411      +    no    1 0.43    0 165.5 1.5e-26   !           -
#> 16      257    410      +    no    1 0.43    0 165.5 1.5e-26   !           -
#> 17      256    409      +    no    1 0.43    0 165.5 1.5e-26   !           -
#> 18      257    410      +    no    1 0.43    0 165.5 1.5e-26   !           -
#> 19      157    310      +    no    1 0.46    0 165.3 1.6e-26   !           -
#> 20      154    307      +    no    1 0.47    0 164.7 2.0e-26   !           -
#> 21      256    409      +    no    1 0.42    0 164.7 2.1e-26   !           -
#> 22      256    409      +    no    1 0.44    0 164.1 2.5e-26   !           -
#> 23      257    410      +    no    1 0.44    0 163.9 2.7e-26   !           -
#> 24      256    409      +    no    1 0.44    0 163.9 2.7e-26   !           -
#> 25      257    410      +    no    1 0.44    0 163.7 2.8e-26   !           -
#> 26      255    408      +    no    1 0.44    0 163.6 3.0e-26   !           -
#> 27      216    368      +    no    1 0.46    0 162.5 4.5e-26   !           -
#> 28      178    331      +    no    1 0.47    0 162.4 4.6e-26   !           -
#> 29      258    411      +    no    1 0.42    0 161.0 7.7e-26   !           -
#> 30      156    309      +    no    1 0.47    0 160.6 8.7e-26   !           -
#> 31      195    348      +    no    1 0.48    0 160.1 1.0e-25   !           -
#> 32      188    340      +    no    1 0.50    0 158.7 1.7e-25   !           -
#> 33      223    376      +    no    1 0.50    0 156.6 3.7e-25   !           -
#> 34      230    383      +    no    1 0.49    0 156.6 3.7e-25   !           -
#> 35      219    371      +    no    1 0.45    0 156.5 3.8e-25   !           -
#> 36      263    416      +    no    1 0.44    0 155.7 5.0e-25   !           -
#> 37      256    408      +    no    1 0.43    0 154.0 9.2e-25   !           -
#> 38      194    346      +    no    1 0.48    0 153.9 9.6e-25   !           -
#> 39      115    268      +    no    1 0.45    0 151.2 2.5e-24   !           -
#> 40      258    414      +    no    1 0.44    0 143.4 4.1e-23   !           -
#> 41      319    475      +    no    1 0.45    0 139.6 1.6e-22   !           -
#> 42      188    340      +    no    1 0.46    0 121.1 1.2e-19   !           -
#> 43      234    388      +    no    1 0.57    0 119.6 2.1e-19   !           -
#> 44       95    248      +    no    1 0.51    0 118.1 3.6e-19   !           -
#> 45       94    246      +    no    1 0.50    0 115.1 1.0e-18   !           -
#> 46       95    246      +    no    1 0.50    0 108.2 1.3e-17   !           -
#> 47       93    246      +    no    1 0.44    0 102.4 9.9e-17   !           -
#> 48      218    373      +    no    1 0.51    0  93.5 2.4e-15   !           -
#> 49        1     10      +    5'    2 0.50    0  -5.4 5.9e+00   ?           -
```

`cmsearch`, by default, returns a table with information about each hit.
However, it can optionally also output an alignment of the hits in
Stockholm format.

``` r
alnfile <- tempfile("alignment-", fileext = ".stk")
cmsearch(cm = cm, seq = sampseqs, alignment = alnfile)
#>    target_name target_accession query_name query_accession mdl mdl_from mdl_to
#> 1        seq45                -  5_8S_rRNA         RF00002  cm        1    154
#> 2         seq3                -  5_8S_rRNA         RF00002  cm        1    154
#> 3         seq2                -  5_8S_rRNA         RF00002  cm        1    154
#> 4        seq28                -  5_8S_rRNA         RF00002  cm        1    154
#> 5        seq23                -  5_8S_rRNA         RF00002  cm        1    154
#> 6         seq9                -  5_8S_rRNA         RF00002  cm        1    154
#> 7         seq7                -  5_8S_rRNA         RF00002  cm        1    154
#> 8        seq48                -  5_8S_rRNA         RF00002  cm        1    154
#> 9         seq5                -  5_8S_rRNA         RF00002  cm        1    154
#> 10        seq6                -  5_8S_rRNA         RF00002  cm        1    154
#> 11       seq11                -  5_8S_rRNA         RF00002  cm        1    154
#> 12       seq12                -  5_8S_rRNA         RF00002  cm        1    154
#> 13       seq13                -  5_8S_rRNA         RF00002  cm        1    154
#> 14       seq17                -  5_8S_rRNA         RF00002  cm        1    154
#> 15       seq18                -  5_8S_rRNA         RF00002  cm        1    154
#> 16       seq19                -  5_8S_rRNA         RF00002  cm        1    154
#> 17       seq25                -  5_8S_rRNA         RF00002  cm        1    154
#> 18       seq47                -  5_8S_rRNA         RF00002  cm        1    154
#> 19       seq14                -  5_8S_rRNA         RF00002  cm        1    154
#> 20        seq8                -  5_8S_rRNA         RF00002  cm        1    154
#> 21       seq21                -  5_8S_rRNA         RF00002  cm        1    154
#> 22       seq36                -  5_8S_rRNA         RF00002  cm        1    154
#> 23       seq10                -  5_8S_rRNA         RF00002  cm        1    154
#> 24       seq22                -  5_8S_rRNA         RF00002  cm        1    154
#> 25       seq50                -  5_8S_rRNA         RF00002  cm        1    154
#> 26       seq44                -  5_8S_rRNA         RF00002  cm        1    154
#> 27       seq27                -  5_8S_rRNA         RF00002  cm        1    154
#> 28       seq35                -  5_8S_rRNA         RF00002  cm        1    154
#> 29       seq15                -  5_8S_rRNA         RF00002  cm        1    154
#> 30       seq38                -  5_8S_rRNA         RF00002  cm        1    154
#> 31       seq30                -  5_8S_rRNA         RF00002  cm        1    154
#> 32       seq26                -  5_8S_rRNA         RF00002  cm        1    154
#> 33       seq33                -  5_8S_rRNA         RF00002  cm        1    154
#> 34       seq46                -  5_8S_rRNA         RF00002  cm        1    154
#> 35       seq34                -  5_8S_rRNA         RF00002  cm        1    154
#> 36       seq39                -  5_8S_rRNA         RF00002  cm        1    154
#> 37       seq16                -  5_8S_rRNA         RF00002  cm        1    154
#> 38       seq31                -  5_8S_rRNA         RF00002  cm        1    154
#> 39        seq1                -  5_8S_rRNA         RF00002  cm        1    154
#> 40       seq41                -  5_8S_rRNA         RF00002  cm        1    154
#> 41       seq29                -  5_8S_rRNA         RF00002  cm        1    154
#> 42       seq49                -  5_8S_rRNA         RF00002  cm        1    154
#> 43       seq20                -  5_8S_rRNA         RF00002  cm        1    154
#> 44        seq4                -  5_8S_rRNA         RF00002  cm        1    154
#> 45       seq32                -  5_8S_rRNA         RF00002  cm        1    154
#> 46       seq24                -  5_8S_rRNA         RF00002  cm        1    154
#> 47       seq37                -  5_8S_rRNA         RF00002  cm        1    154
#> 48       seq43                -  5_8S_rRNA         RF00002  cm        1    154
#> 49       seq37                -  5_8S_rRNA         RF00002  cm      145    154
#>    seq_from seq_to strand trunc pass   gc bias score E_value inc description
#> 1       295    448      +    no    1 0.47    0 171.3 1.9e-27   !           -
#> 2       236    389      +    no    1 0.45    0 171.0 2.1e-27   !           -
#> 3       193    346      +    no    1 0.48    0 170.7 2.4e-27   !           -
#> 4       192    345      +    no    1 0.48    0 170.7 2.4e-27   !           -
#> 5       194    347      +    no    1 0.47    0 168.5 5.3e-27   !           -
#> 6       170    323      +    no    1 0.47    0 167.3 7.9e-27   !           -
#> 7       191    344      +    no    1 0.50    0 167.3 8.1e-27   !           -
#> 8       192    345      +    no    1 0.49    0 167.1 8.7e-27   !           -
#> 9       257    410      +    no    1 0.43    0 165.5 1.5e-26   !           -
#> 10      257    410      +    no    1 0.43    0 165.5 1.5e-26   !           -
#> 11      258    411      +    no    1 0.43    0 165.5 1.5e-26   !           -
#> 12      257    410      +    no    1 0.43    0 165.5 1.5e-26   !           -
#> 13      256    409      +    no    1 0.43    0 165.5 1.5e-26   !           -
#> 14      257    410      +    no    1 0.43    0 165.5 1.5e-26   !           -
#> 15      258    411      +    no    1 0.43    0 165.5 1.5e-26   !           -
#> 16      257    410      +    no    1 0.43    0 165.5 1.5e-26   !           -
#> 17      256    409      +    no    1 0.43    0 165.5 1.5e-26   !           -
#> 18      257    410      +    no    1 0.43    0 165.5 1.5e-26   !           -
#> 19      157    310      +    no    1 0.46    0 165.3 1.6e-26   !           -
#> 20      154    307      +    no    1 0.47    0 164.7 2.0e-26   !           -
#> 21      256    409      +    no    1 0.42    0 164.7 2.1e-26   !           -
#> 22      256    409      +    no    1 0.44    0 164.1 2.5e-26   !           -
#> 23      257    410      +    no    1 0.44    0 163.9 2.7e-26   !           -
#> 24      256    409      +    no    1 0.44    0 163.9 2.7e-26   !           -
#> 25      257    410      +    no    1 0.44    0 163.7 2.8e-26   !           -
#> 26      255    408      +    no    1 0.44    0 163.6 3.0e-26   !           -
#> 27      216    368      +    no    1 0.46    0 162.5 4.5e-26   !           -
#> 28      178    331      +    no    1 0.47    0 162.4 4.6e-26   !           -
#> 29      258    411      +    no    1 0.42    0 161.0 7.7e-26   !           -
#> 30      156    309      +    no    1 0.47    0 160.6 8.7e-26   !           -
#> 31      195    348      +    no    1 0.48    0 160.1 1.0e-25   !           -
#> 32      188    340      +    no    1 0.50    0 158.7 1.7e-25   !           -
#> 33      223    376      +    no    1 0.50    0 156.6 3.7e-25   !           -
#> 34      230    383      +    no    1 0.49    0 156.6 3.7e-25   !           -
#> 35      219    371      +    no    1 0.45    0 156.5 3.8e-25   !           -
#> 36      263    416      +    no    1 0.44    0 155.7 5.0e-25   !           -
#> 37      256    408      +    no    1 0.43    0 154.0 9.2e-25   !           -
#> 38      194    346      +    no    1 0.48    0 153.9 9.6e-25   !           -
#> 39      115    268      +    no    1 0.45    0 151.2 2.5e-24   !           -
#> 40      258    414      +    no    1 0.44    0 143.4 4.1e-23   !           -
#> 41      319    475      +    no    1 0.45    0 139.6 1.6e-22   !           -
#> 42      188    340      +    no    1 0.46    0 121.1 1.2e-19   !           -
#> 43      234    388      +    no    1 0.57    0 119.6 2.1e-19   !           -
#> 44       95    248      +    no    1 0.51    0 118.1 3.6e-19   !           -
#> 45       94    246      +    no    1 0.50    0 115.1 1.0e-18   !           -
#> 46       95    246      +    no    1 0.50    0 108.2 1.3e-17   !           -
#> 47       93    246      +    no    1 0.44    0 102.4 9.9e-17   !           -
#> 48      218    373      +    no    1 0.51    0  93.5 2.4e-15   !           -
#> 49        1     10      +    5'    2 0.50    0  -5.4 5.9e+00   ?           -
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
defined by the CM. These are given as column (“GC”) annotations in the
Stockholm alignment file:

``` r
msa$GC$SS_cons
#> 164-letter BString object
#> seq: ::::::::::::::::::::::::::::::::::::...<<..____.>>>>>>>>>::::::::::::::::::
msa$GC$RF
#> 164-letter BString object
#> seq: AACuuUuAgCGAUGGAUguCUuGGCUCccGuaUCGA...gg..Uuuu.cccgggggCAUgccUGuuugAGUGUCa
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
#> RNAStringSet object of length 48:
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
aln$GC$SS_cons
#> 164-letter BString object
#> seq: ::::::::::::::::::::::::::::::::::::...<<..____.>>>>>>>>>::::::::::::::::::
aln$GC$RF
#> 164-letter BString object
#> seq: AACuuUuAgCGAUGGAUguCUuGGCUCccGuaUCGA...gg..Uuuu.cccgggggCAUgccUGuuugAGUGUCa
```

## cmbuild

`cmbuild` is used to create new CMs from annotated multiple sequence
alignments. To illustrate the process, we use the seed alignment for the
5.8S rRNA CM from RFAM. It is included as a sample file in `inferrnal`.

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
