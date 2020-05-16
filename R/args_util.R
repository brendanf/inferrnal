# pipe-friendly utility functions to add command line arguments of
# various types.

arg2opt <- function(arg) {
    assertthat::assert_that(assertthat::is.string(arg), nchar(arg) > 0)
    paste0(ifelse(nchar(arg) == 1, "-", "--"), arg)
}

check_null_args <- function(argval, argname, allow_null) {
    if (is.null(argval) || is.na(argval)) {
        if (isTRUE(allow_null)) {
            return(TRUE)
        } else {
            assertthat::assert_that(
                assertthat::is.flag(allow_null),
                !is.na(allow_null)
            )
            stop(
                if (is.null(argval)) "NULL" else "NA",
                " value for option '", arg2opt(argname),
                "' but allow_null is FALSE."
            )
        }
    }
    FALSE
}

flag_opt <- function(argval, argname = deparse(substitute(argval)),
                    invert = FALSE) {
    if (isTRUE(argval)) {
        if (isFALSE(invert)) {
            return(c(arg2opt(argname)))
        } else if (isTRUE(invert)) {
            return(character())
        }
    } else if (isFALSE(argval)) {
        if (isTRUE(invert)) {
            return(c(arg2opt(argname)))
        } else if (isFALSE(invert)) {
            return(character())
        }
    }
    # this should always fail, but it will explain why in a standard way.
    assertthat::assert_that(
        assertthat::is.flag(argval),
        !is.na(argval),
        assertthat::is.flag(argval),
        !is.na(invert)
    )
}

multiflag_opt <- function(argval, choices, allow_null = TRUE) {
    assertthat::assert_that(
        is.character(choices),
        assertthat::is.flag(allow_null),
        !is.na(allow_null)
    )
    if (is.null(argval) || is.na(argval)) {
        if (isTRUE(allow_null)) {
            return(character())
        } else {
            stop(
                if (is.null(argval)) "NULL" else "NA",
                " value for multiflag option but 'allow_null' is FALSE."
            )
        }
    }
    assertthat::assert_that(
        assertthat::is.string(argval),
        argval %in% choices
    )
    c(arg2opt(argval))
}

string_opt <- function(argval, argname = deparse(substitute(argval)),
                        choices = NULL, allow_null = TRUE) {
    if (check_null_args(argval, argname, allow_null)) return(character())
    assertthat::assert_that(
        assertthat::is.string(argval)
    )
    if (!is.null(choices)) {
        assertthat::assert_that(
            is.character(choices),
            argval %in% choices
        )
    }
    c(arg2opt(argname), argval)
}

infile_opt <- function(argval, argname = deparse(substitute(argval)),
                        allow_null = TRUE) {
    assertthat::assert_that(assertthat::is.string(argname))
    if (check_null_args(argval, argname, allow_null)) return(character())
    assertthat::assert_that(
        assertthat::is.string(argval),
        assertthat::is.readable(argval)
    )
    c(arg2opt(argname), argval)
}

float_opt <- function(argval, argname = deparse(substitute(argval)),
                        allow_null = TRUE) {
    assertthat::assert_that(assertthat::is.string(argname))
    if (check_null_args(argval, argname, allow_null)) return(character())
    assertthat::assert_that(assertthat::is.number(argval))
    c(arg2opt(argname), format(argval, scientific = FALSE))
}

nonneg_float_opt <- function(argval, argname = deparse(substitute(argval)),
                            allow_null = TRUE) {
    assertthat::assert_that(assertthat::is.string(argname))
    if (check_null_args(argval, argname, allow_null)) return(character())
    assertthat::assert_that(
        assertthat::is.number(argval),
        argval >= 0
    )
    c(arg2opt(argname), format(argval, scientific = FALSE))
}

fraction_opt <- function(argval, argname = deparse(substitute(argval)),
                        allow_null = TRUE) {
    assertthat::assert_that(assertthat::is.string(argname))
    if (check_null_args(argval, argname, allow_null)) return(character())
    assertthat::assert_that(
        assertthat::is.number(argval),
        argval >= 0,
        argval <= 1
    )
    c(arg2opt(argname), format(argval, scientific = FALSE))
}

percent_opt <- function(argval, argname = deparse(substitute(argval)),
                        allow_null = TRUE) {
    assertthat::assert_that(assertthat::is.string(argname))
    if (check_null_args(argval, argname, allow_null)) return(character())
    assertthat::assert_that(
        assertthat::is.number(argval),
        argval >= 0,
        argval <= 100
    )
    c(arg2opt(argname), format(argval, scientific = FALSE))
}

integer_opt <- function(argval, argname = deparse(substitute(argval)),
                        allow_null = TRUE) {
    assertthat::assert_that(assertthat::is.string(argname))
    if (check_null_args(argval, argname, allow_null)) return(character())
    assertthat::assert_that(
        assertthat::is.number(argval),
        abs(argval - round(argval)) < .Machine$double.eps^0.5
        )
    c(arg2opt(argname), format(round(argval), scientific = FALSE))
}

count_opt <- function(argval, argname = deparse(substitute(argval)),
                        allow_null = TRUE) {
    assertthat::assert_that(assertthat::is.string(argname))
    if (check_null_args(argval, argname, allow_null)) return(character())
    assertthat::assert_that(assertthat::is.count(argval))
    c(arg2opt(argname), format(round(argval), scientific = FALSE))
}
