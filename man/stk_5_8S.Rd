% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/example_data.R
\name{stk_5_8S}
\alias{stk_5_8S}
\title{Annotated Seed Alignment for 5.8S rRNA}
\usage{
stk_5_8S()
}
\value{
(\code{character} string) the path to the stk file
}
\description{
This is a convenience function to return the path the seed alignment for the
RFAM 5.8S rRNA covariance model,
for use in \code{\link[=cmbuild]{cmbuild()}}.
The original file is from
\url{https://rfam.xfam.org/family/RF00002/stockholm}.
}
\examples{
    cmbuild(
        msafile = stk_5_8S(),
        cmfile_out = "/dev/null",
        force = TRUE,
        quiet = FALSE
    )
}
