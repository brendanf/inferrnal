# Changes in version 0.99.5

* Correct argument `"glocal"` of `cmalign` to `"global"`.  Deprecate `"glocal"`.
* Implement `cmbuild()`, along with all basic options.
* Implement all basic options for `cmsearch()` and `cmalign()`.
* Add an `"extra"` argument to all three `cm*()` functions to allow passing
  advanced options.
* Add a `"quiet"` argument to suppress the sometimes very large console output
  of `cmsearch()`.
* `read_stockholm_msa()` parses and returns residue (GR), sequence (GS), and
  file (GF) annotations.
* `cmalign()` returns an error rather than an empty alignment when Infernal
  exits with an error.
* If `read_stockholm_msa()` is given an open connection, it does not close it
  when done. (But if given a filename or an unopened connection, it opens and
  then closes.)

# Changes in version 0.99.4

* Properly handle explicit `NULL` values for `cpu` and `mxsize`.

# Changes in version 0.99.3

* Added `mxsize` argument to `cmalign` to change the maximum matrix size.

# Changes in version 0.99.2

* Added convenience functions `cm_5_8S()`, `sample_rRNA_fasta()`,
  `sample_rRNA_5_8S()`, and `sample_rRNA_stk()` to access example data files.

# Changes in version 0.99.1

* Added a `NEWS.md` file to track changes to the package.
* Added support for objects of class `ShortRead` as input.
