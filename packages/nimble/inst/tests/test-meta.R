context("Basic checks to make sure the testing itself is working")

a <- 5

test_that("A test to see the testing system is working. This test should pass.",
         expect_that(a*5, equals(25)))

## removed (intentionally) failing test
##test_that("A test to see the testing system is working. This test should fail.",
##         expect_that(a*5, equals(26)))


