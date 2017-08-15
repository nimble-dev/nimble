# This file consists of unit tests of the code generation internals of NIMBLE, and
# should never invoke a compiler.

source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Unit tests of code transformations")

# These wrappers help convert between exprClass and R parse trees.
code_to_expr <- nimble:::RparseTree2ExprClasses
expr_to_code <- function(expr) {
    if (!inherits(expr, 'exprClass')) stop(paste('Expected an exprClass, got', class(expr)))
    deparsed <- nimble:::nimDeparse(expr)
    # nimDeparse returns a list rather than a character vector,
    # and after parsing, the code is wrapped in expression(-).
    code <- parse(text = as.character(deparsed))[[1]]
    return(code)
}

test_that('Transform from code to expr and back', {
    expected_code <- quote({
        x <- rep(0, 4)
        y <- sin(inprod(x, x))
    })
    expr <- code_to_expr(expected_code)
    actual_code <- expr_to_code(expr)
    expect_equal(actual_code, expected_code)
})

test_that('sin(inprod(x,x)) is transformed to sin(eigenInprod(x,x))', {
    code <- quote({
        x <- rep(0, 4)
        y <- sin(inprod(x, x))
    })
    # TODO(fritzo,perrydv) Refactor to expose compiler stages, e.g.
    # - Try to avoid accessing mutable state.
    # - Prefer pure functions over class methods with side-effects.
    # - Add helpers to convert and compare exprClass instances.
    skip('This does not work yet')
    proc <- RCfunProcessing$new()
    compileInfo <- nimble:::RCfunctionCompileClass$new(origRcode = code, newRcode = code)
    expr <- code_to_expr(code)
    expect_equal(actual_code, code)
})
