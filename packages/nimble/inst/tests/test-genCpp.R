# This file consists of unit tests of the code generation internals of NIMBLE, and
# should never invoke a compiler.

source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

skip_if_not_installed('pryr')

context("Unit tests of code transformations")

RparseTree2ExprClasses <- nimble:::RparseTree2ExprClasses
nimDeparse <- nimble:::nimDeparse

test_that('nimDeparse roughly inverts RparseTree2ExprClasses', {
    code <- quote({
        x <- rep(0, 4)
        y <- sin(inprod(x, x))
    })
    nimCode <- RparseTree2ExprClasses(code)
    deparsed <- nimDeparse(nimCode)
    # nimDeparse returns a list rather than a character vector,
    # and the code is wrapped in expression(-).
    actual_code <- parse(text = as.character(deparsed))[[1]]
    expect_equal(actual_code, code)
})
