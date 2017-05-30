# Guidelines for Contributing

## Opening an issue

To file a bug or request a feature,
[open an issue on Github](https://github.com/nimble-dev/nimble/issues/new).

## Submitting code via pull requests

To submit code:

1.  Create a branch. Branches should be named like `my-descriptive-branch`.
2.  Publish your branch to the NIMBLE repo on Github.
3.  Create a pull request. Write a clear description; request review; and
    label the PR as *enhancement*, *bug*, or *cleanup*.
4.  When a reviewer comments **"LGTM"** (looks good to me) and
    Travis CI tests pass, you can merge your PR.
    Prefer **Squash and Merge** over **Merge** unless your git history is
    exceptionally clean and each commit passes CI tests.

## Testing

### Test all features and bugfixes

Add tests for every new feature.
If your new feature is large, consider adding a new test file named `$REPO/packages/nimble/inst/tests/test-<my-new-feature>.R`.

Add a regression test for every bug you fix.

### Wrap tests in `test_that()`

We use Hadley Wickham's [`testthat`](https://www.rdocumentation.org/packages/testthat/versions/1.0.2) package.

Wrap each test in a [`test_that()`](https://www.rdocumentation.org/packages/testthat/versions/1.0.2/topics/test_that) that looks like
```diff
+ # GOOD Wraps a test in test_that().
  test_that('my new feature behaves correctly', {
      # ...test code here...
  })
```
This pattern has many advantages.
First, it allows all tests in a `test-*.R` file to be run, even if one of the tests fails.
The `testthat` package will run all tests and then report all failing tests.
Second, wrapping your test in a `test_that(_, {_})` prevents tests from leaking variables into the global namespace.
This is especially important when you copy-and-modify to create tests, for example
```diff
- # BAD Global definitions allows uncaught errors.
  fun42 <- nimbleFunction(run = function(x = double(0)) {
      returnType(double(0))
      return(x + x)
  })
  expect_equal(fun42(0), compileNimble(fun42)(0))
  
  fun43 <- nimbleFunction(run = function(x = double(0)) {
      returnType(double(0))
      return(x * x)
  })
  expect_equal(fun43(0), compileNimble(fun42)(0))  # <-- Uncaught typo here.

+ # GOOD Local definitions avoid possibility of errors.
  test_that('Addition works', {
      fun <- nimbleFunction(run = function(x = double(0)) {
          returnType(double(0))
          return(x + x)
      })
      expect_equal(fun(0), compileNimble(fun)(0))
  })
  
  test_that('Multiplication works', {
      fun <- nimbleFunction(run = function(x = double(0)) {
          returnType(double(0))
          return(x * x)
      })
      expect_equal(fun(0), compileNimble(fun)(0))   # <-- No possibility of typo.
  })
```
Note that `compileNimble()` requires all user-defined functions to exist in the global environment; use the [`temporarilyAssignInGlobalEnv()`](https://github.com/nimble-dev/nimble/blob/b4129f2/packages/nimble/inst/tests/test_utils.R#L22) function to do this cleanly from within a `test_that()` block, e.g.
```diff
+ # GOOD Use temporarilyAssignInGlobalEnv() inside test_that().
  test_that('One RCfunction can use another', {
      helper <- nimbleFunction(run = function(){
          returnType(double(0))
          return(3.1415)
      })
      temporarilyAssignInGlobalNamespace(helper)  # <-- Make helper available.
      fun <- nimbleFunction(run = function(r = double(0)){
          returnType(double(0))
          return(helper() * r * r)
      })
      compileNimble(fun)
  })
```

### Use parameterized tests

Use parametrized tests to test many possibilities in a for loop.
We recommended two different practices for parametrized testing: `for`-outside-`test_that` and `for`-inside-`test_that`.
```diff
+ # GOOD Create multiple test_that()s inside a for loop.
  sizes <- c(1,2,3,10,100)
  for (size in sizes) {
      test_that(paste('Function works with vector of size', size), {
          x <- rep(1, size)
          y <- rep(2, size)
          expect_equal(x + x, y)
      })
  }

+ # GOOD Embed multiple expect()s inside a single test_that()
+ # and show parameter values on error using the info arg of expect_equal().
  test_that('Function works with vectors', {
      sizes <- c(1,2,3,10,100)
      for (size in sizes) {
          x <- rep(1, size)
          y <- rep(2, size)
          expect_equal(x + x, y, info = paste(' where size =', size))
      })
  }
  
- # BAD Embed multiple expects() in a for loop, but forget to note parameter:
- # errors will fail to print which `size` caused the error.
  test_that('Function works with vectors', {
      sizes <- c(1,2,3,10,100)
      for (size in sizes) {
          x <- rep(1, size)
          y <- rep(2, size)
          expect_equal(x + x, y)  # <-- Missing info arg.
      })
  }
```

### Test expected errors with `expect_error`

Use `expect_error` to test for expected errors,
e.g. when testing that our compiler correctly errors on some code.
```diff
+ # GOOD use expect_error
  test_that('Compiler fails when no returnType is provided', {
      nf <- nimbleFunction(run = function(){ return(0) })
      expect_error(compileNimble(nf))
  })
```

### Use `test_utils.R`

See also the numerous testing tools in [`test_utils.R`](https://github.com/nimble-dev/nimble/blob/devel/packages/nimble/inst/tests/test_utils.R), including:

-   `system.in.dir` -
    Runs a command in a directory in a platform independent way.
-   `temporarilyAssignInGlobalEnv` -
    Temporarily assigns a variable in the global environment, as needed by
    user-defined functions that are needed by other nimble functions.
-   `test_mcmc` - A function for testing MCMC (called from `test_mcmc.R`
-   `clearOldOutput` - Remove a file if it does not exist (like `rm -f`).

and many undocumented test utilities.

### Run tests with `test_package`

Run tests locally using e.g. `test_package('nimble', 'my-component')` to test `$REPO/packages/nimble/inst/tests/test-my-component.R`.
Note that although can be using either `test_package()` or directly via `Rscript`, the `Rscript` method fails to detect some errors.

## Style

> Programs must be written for people to read,
> and only incidentally for machines to execute.
>
> --<cite>Harold Abelson, Structure and Interpretation of Computer Programs</cite>

Try to be consistent with the code you're modifying.
If in doubt, seek guidance from
[Google's R Style Guide](https://google.github.io/styleguide/Rguide.xml) or
[Hadley Wickham's style guide](http://adv-r.had.co.nz/Style.html).

### Write clear comments

Write clear concise comments aimed at other developers.
If you write a great internal helper function,
and you want other developers to be able to use your code,
then write comment describing what your helper is and how to use it.
Write comments grammatically:
sentences should start with a capital letter and end with a period.

### Avoid commented-out code

Avoid submitting commented-out code into the repo.
If you really want to remember the code, then remove it from the repo with a descriptive commit message.
If you really want the code to exist, then gate it by an `if()else` with an argument or global option.
