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


test_that("can read file", {
    expect_known_value(
        read_stockholm_msa(
            system.file(
                file.path("extdata", "sample.stk"),
                package = "inferrnal"
            )
        ),
        file = "read_stockholm_msa_out"
    )
})
