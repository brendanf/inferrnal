test_that("arg2opt works", {
    expect_error(arg2opt(3))
    expect_error(arg2opt(NULL))
    expect_error(arg2opt(NA_character_))
    expect_error(arg2opt(""))
    expect_identical(arg2opt("a"), "-a")
    expect_identical(arg2opt("option"), "--option")
})

test_that("flag_opt works", {
    expect_error(flag_opt("true", "a"))
    expect_error(flag_opt(NA, "a"))
    expect_error(flag_opt(TRUE, "a", invert = NA))
    expect_identical(
        flag_opt(TRUE, "a"),
        c("-a")
    )
    expect_identical(
        flag_opt(FALSE, "a"),
        character()
    )
    expect_identical(
        flag_opt(FALSE, "a", invert = TRUE),
        c("-a")
    )
    expect_identical(
        flag_opt(TRUE, "a", invert = TRUE),
        character()
    )
})

test_that("multiflag_opt works", {
    expect_error(multiflag_opt("true", c("a", "b", "c")))
    expect_error(multiflag_opt(NA, c("a", "b", "c"), allow_null = FALSE))
    expect_error(multiflag_opt(TRUE, c("a", "b", "c")))
    expect_error(multiflag_opt(10, c("a", "b", "c")))
    expect_error(multiflag_opt(10L, c("a", "b", "c")))
    expect_error(multiflag_opt(NULL, c("a", "b", "c"), allow_null = FALSE))
    expect_error(multiflag_opt("b", "a"))
    expect_identical(
        multiflag_opt(NULL, c("a", "b", "c"), allow_null = TRUE),
        character()
    )
    expect_identical(
        multiflag_opt("a", c("a", "b", "c")),
        c("-a")
    )
    expect_identical(
        multiflag_opt("b", c("a", "b", "c")),
        c("-b")
    )
    expect_identical(
        multiflag_opt("a", c("a", "b", "c"), allow_null = TRUE),
        c("-a")
    )
})

test_that("string_opt works", {
    expect_error(string_opt("d",    "a", c("a", "b", "c")))
    expect_error(string_opt(NA,     "a", c("a", "b", "c"), allow_null = FALSE))
    expect_error(string_opt(TRUE,   "a", c("a", "b", "c")))
    expect_error(string_opt(10,     "a", c("a", "b", "c")))
    expect_error(string_opt(10L,    "a", c("a", "b", "c")))
    expect_error(string_opt(NULL,   "a", c("a", "b", "c"), allow_null = FALSE))
    expect_identical(
        string_opt(NULL, "a", allow_null = TRUE),
        character()
    )
    expect_identical(
        string_opt(NULL, "a", c("a", "b", "c"), allow_null = TRUE),
        character()
    )
    expect_identical(
        string_opt("a", "d", c("a", "b", "c")),
        c("-d", "a")
    )
    expect_identical(
        string_opt("a", "d"),
        c("-d", "a")
    )
})

temp <- tempfile()
writeLines("test", temp)

test_that("infile_opt works", {
    expect_error(infile_opt("d",    "a"))
    expect_error(infile_opt(TRUE,   "a"))
    expect_error(infile_opt(10,     "a"))
    expect_error(infile_opt(10L,    "a"))
    expect_error(infile_opt(NA,     "a", allow_null = FALSE))
    expect_error(infile_opt(NULL,   "a", allow_null = FALSE))
    expect_identical(
        string_opt(temp, "a", allow_null = TRUE),
        c("-a", temp)
    )
    expect_identical(
        string_opt(NULL, "a", allow_null = TRUE),
        character()
    )
})

test_that("float_opt works", {
    expect_error(float_opt("d",  "a"))
    expect_error(float_opt(TRUE, "a"))
    expect_error(float_opt(NULL, "a", allow_null = FALSE))
    expect_error(float_opt(NA,   "a", allow_null = FALSE))
    expect_identical(
        float_opt(NULL, "a", allow_null = TRUE),
        character()
    )
    expect_identical(
        float_opt(NA, "a", allow_null = TRUE),
        character()
    )
    expect_identical(
        float_opt(0.000000001, "d", allow_null = TRUE),
        c("-d", "0.000000001")
    )
    expect_identical(
        float_opt(0.000000001, "d"),
        c("-d", "0.000000001")
    )
})

test_that("nonneg_float_opt works", {
    expect_error(nonneg_float_opt("d",  "a"))
    expect_error(nonneg_float_opt(TRUE, "a"))
    expect_error(nonneg_float_opt(NULL, "a", allow_null = FALSE))
    expect_error(nonneg_float_opt(NA,   "a", allow_null = FALSE))
    expect_error(nonneg_float_opt(-1,   "a", allow_null = FALSE))
    expect_error(nonneg_float_opt(-1,   "a"))
    expect_identical(
        nonneg_float_opt(NULL, "a", allow_null = TRUE),
        character()
    )
    expect_identical(
        nonneg_float_opt(NA, "a", allow_null = TRUE),
        character()
    )
    expect_identical(
        nonneg_float_opt(0.000000001, "d", allow_null = TRUE),
        c("-d", "0.000000001")
    )
    expect_identical(
        nonneg_float_opt(0.000000001, "d"),
        c("-d", "0.000000001")
    )
})

test_that("fraction_opt works", {
    expect_error(fraction_opt("d",  "a"))
    expect_error(fraction_opt(TRUE, "a"))
    expect_error(fraction_opt(NULL, "a", allow_null = FALSE))
    expect_error(fraction_opt(NA,   "a", allow_null = FALSE))
    expect_error(fraction_opt(-1,   "a", allow_null = FALSE))
    expect_error(fraction_opt(-1,   "a"))
    expect_error(fraction_opt(1.01,   "a", allow_null = FALSE))
    expect_error(fraction_opt(1.01,   "a"))
    expect_identical(
        fraction_opt(NULL, "a", allow_null = TRUE),
        character()
    )
    expect_identical(
        fraction_opt(NA, "a", allow_null = TRUE),
        character()
    )
    expect_identical(
        fraction_opt(0.000000001, "d", allow_null = TRUE),
        c("-d", "0.000000001")
    )
    expect_identical(
        fraction_opt(0.000000001, "d"),
        c("-d", "0.000000001")
    )
})

test_that("percent_opt works", {
    expect_error(percent_opt("d",  "a"))
    expect_error(percent_opt(TRUE, "a"))
    expect_error(percent_opt(NULL, "a", allow_null = FALSE))
    expect_error(percent_opt(NA,   "a", allow_null = FALSE))
    expect_error(percent_opt(-1,   "a", allow_null = FALSE))
    expect_error(percent_opt(-1,   "a"))
    expect_error(percent_opt(101,   "a", allow_null = FALSE))
    expect_error(percent_opt(101,   "a"))
    expect_identical(
        percent_opt(NULL, "a", allow_null = TRUE),
        character()
    )
    expect_identical(
        percent_opt(NA, "a", allow_null = TRUE),
        character()
    )
    expect_identical(
        percent_opt(0.000000001, "d", allow_null = TRUE),
        c("-d", "0.000000001")
    )
    expect_identical(
        percent_opt(0.000000001, "d"),
        c("-d", "0.000000001")
    )
})

test_that("integer_opt works", {
    expect_error(integer_opt("d",  "a"))
    expect_error(integer_opt(TRUE, "a"))
    expect_error(integer_opt(NULL, "a", allow_null = FALSE))
    expect_error(integer_opt(NA,   "a", allow_null = FALSE))
    expect_error(integer_opt(1.5,   "a"))
    expect_identical(
        integer_opt(NULL, "a", allow_null = TRUE),
        character()
    )
    expect_identical(
        integer_opt(NA, "a", allow_null = TRUE),
        character()
    )
    expect_identical(
        integer_opt(1, "d", allow_null = TRUE),
        c("-d", "1")
    )
    expect_identical(
        integer_opt(-1, "d"),
        c("-d", "-1")
    )
    expect_identical(
        integer_opt(1L, "d", allow_null = TRUE),
        c("-d", "1")
    )
    expect_identical(
        integer_opt(-1L, "d"),
        c("-d", "-1")
    )
})



test_that("count_opt works", {
    expect_error(count_opt("d",  "a"))
    expect_error(count_opt(TRUE, "a"))
    expect_error(count_opt(NULL, "a", allow_null = FALSE))
    expect_error(count_opt(NA,   "a", allow_null = FALSE))
    expect_error(count_opt(1.5,   "a"))
    expect_error(count_opt(-1,   "a"))
    expect_error(count_opt(-1L,   "a"))
    expect_identical(
        count_opt(NULL, "a", allow_null = TRUE),
        character()
    )
    expect_identical(
        count_opt(NA, "a", allow_null = TRUE),
        character()
    )
    expect_identical(
        count_opt(1, "d", allow_null = TRUE),
        c("-d", "1")
    )
    expect_identical(
        count_opt(1, "d"),
        c("-d", "1")
    )
})

b = TRUE
test_flag = TRUE
s = "b1"
test_str = "b1"
f = temp
test_file = temp
n = 10
test_int = 10
x = 0.5
test_num = 0.5
test_that("NSE for option names works", {
    expect_identical(flag_opt(b), c("-b"))
    expect_identical(flag_opt(test_flag), c("--test_flag"))
    expect_identical(string_opt(s), c("-s", s))
    expect_identical(string_opt(test_str), c("--test_str", test_str))
    expect_identical(infile_opt(f), c("-f", f))
    expect_identical(infile_opt(test_file), c("--test_file", test_file))
    expect_identical(float_opt(x), c("-x", x))
    expect_identical(float_opt(test_num), c("--test_num", test_num))
    expect_identical(nonneg_float_opt(x), c("-x", x))
    expect_identical(nonneg_float_opt(test_num), c("--test_num", test_num))
    expect_identical(fraction_opt(x), c("-x", x))
    expect_identical(fraction_opt(test_num), c("--test_num", test_num))
    expect_identical(percent_opt(x), c("-x", x))
    expect_identical(percent_opt(test_num), c("--test_num", test_num))
    expect_identical(integer_opt(n), c("-n", n))
    expect_identical(integer_opt(test_int), c("--test_int", test_int))
    expect_identical(count_opt(n), c("-n", n))
    expect_identical(count_opt(test_int), c("--test_int", test_int))
})
