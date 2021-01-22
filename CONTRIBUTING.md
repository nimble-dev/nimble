# Guidelines for Contributing

Nimble is managed across multiple repos:
[nimble](https://github.com/nimble-dev/nimble) for code and the user manual,
[nimble-docs](https://github.com/nimble-dev/nimble-doc)
for developer documentation,
and miscellaneous other repos under the [nimble-dev](https://github.com/nimble-dev) org.
We welcome issues and pull requests to all of these repos.

## Opening an issue

To file a bug or request a feature or documentation,
[open an issue on Github](https://github.com/nimble-dev/nimble/issues/new).
When possible, include a fully reproducible example demonstrating the error or unwanted behaviour.
Please note your operating system, R version, and NIMBLE version (either CRAN X.X-X or git commit).

```sh
uname -a
R --version
Rscript -e "packageVersion('nimble')"
```

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

## Documentation

To keep documentation in sync with code, try to update the
[User Manual](UserManual/)
in the same pull request as the code that implements a change.
This eases the job of reviewers, and makes it easier for users to use features that are experimental.

To rebuild function help info after updating function documentation, run

```sh
cd packages
make man
```

To rebuild the [NimbleUserManual.pdf](UserManual/NimbleUserManual.pdf)
after adding sections, run

```sh
cd UserManual
make
```

Take care to review the generated `NimbleUserManual.pdf` before submitting, since knitr
won't tell notify about execution errors that happend while building the manual.

### Experimental features

To encourage sharing of code, we try to merge our branches early and often.
This means that some experimental features will be in `devel` and `master`
branches. If you are writing an experimental feature, please document the status
by writing "EXPERIMENTAL" in the feature's documentation.

Often we hide experimental features behind *feature toggles*.
If your feature toggle is temporary or only partially tested, please name it
`experimentalMyFeature` to indicate to users that the feature is still
experimental.

## Testing

We test on travis, but you can get quicker results by running tests locally in parallel.
```sh
cd $REPO/packages
make test           # This runs tests in parallel.
```

### Test all features and bugfixes

Add tests for every new feature.
If your new feature is large, consider adding a new test file named
`$REPO/packages/nimble/tests/testthat/test-<my-new-feature>.R`.

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
Note that `compileNimble()` requires all user-defined functions to exist in the global environment; use the [`temporarilyAssignInGlobalEnv()`](https://github.com/nimble-dev/nimble/blob/b4129f2/packages/nimble/tests/testthat/test_utils.R#L22) function to do this cleanly from within a `test_that()` block, e.g.
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

### Mark known failures with `expect_failure()` or `skip()`

The best way to mark a known failure is to wrap an individually failing expectaion in `expect_failure()`.
The second best way to mark a known failure is to call `skip()` within a `test_that()` block.
**AVOID** the `try(test_that())` idiom, since this behaves differently when running tests under `Rscript` versus `test_package()/test_file()`.
```diff
+ # GOOD Wrap an expectation in expect_failure() and point to a github issue.
  test_that('One equals zero', {
      expect_failure(
          expect_equal(1, 0),
          info = 'KNOWN ISSUE https://github.com/nimble-dev/nimble/issues/10'
      )
  })
  
+ # GOOD skip() all or part of a test and point to a github issue.
  test_that('One equals everything', {
      expect_equal(1, 1)  ## This passes.
      skip('KNOWN ISSUE https://github.com/nimble-dev/nimble/issues/10')
      expect_equal(1, 0)  ## This fails.
  })
  
- # BAD Avoid wrapping tests in try(), since they will still fail under `test_package()`.
- # This passes under Rscript but fails under test_package() and test_file().
  try(test_that('One equals zero', {
      expect_equal(1, 0)
  }))
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

If you are using parametrized tests and some of the tests are known failures of NIMBLE, you can pass in a flag to the parametrized test and query that flag with `wrap_if_matches` so as to only wrap the specific tests that are known to fail with an `expect_failure` (or in some cases `expect_warning`). See `mathTestLists.R` and `test-mcmc.R` for examples. Note that `test-size.R` was developed before the system recommended in this document and that `test-size.R` handles very convoluted cases, so do not follow the syntax used in `test-size.R`.

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

Note that if we would like an error to be triggered but it is not, this would call for wrapping the `expect_error` call in an `expect_failure` call. 

### Use `test_utils.R`

See also the numerous testing tools in [`test_utils.R`](https://github.com/nimble-dev/nimble/blob/devel/packages/nimble/tests/testthat/test_utils.R), including:

-   `system.in.dir` -
    Runs a command in a directory in a platform independent way.
-   `temporarilyAssignInGlobalEnv` -
    Temporarily assigns a variable in the global environment, as needed by
    user-defined functions that are needed by other nimble functions.
-   `test_mcmc` - A function for testing MCMC (called from `test_mcmc.R`
-   `clearOldOutput` - Remove a file if it does not exist (like `rm -f`).

and many undocumented test utilities.

### Using goldfiles and determining what output to include/suppress

As seen in `test-mcmc.R` and others, we compare the printed output (particularly for MCMC summaries) with known results from previous testing to ensure that numeric values have not changed. Please follow the template in `test-mcmc.R`. Goldfiles should be saved in `tests/testthat` and included in the repository and should be updated as needed, but only with care that any changes relative to older goldfiles are understood. 

We recommend including `options(warn=1)` so that warnings are printed out. However, note that it is difficult to ensure that warnings appear in the goldfiles and are able to be checked exactly against previous testing instantiations. In some cases we use `sink_with_messages` rather than `sink` to capture output from testing. This does not capture the `W` values printed out by testthat, but it does (e.g., in the case of `test-models.R`) capture other warning messages that one may want to check versus previous versions of the code.

We recommend turning verbosity and MCMC progress bars off during testing to avoid extensive unimportant output. Again see `test-mcmc.R`. 

### Run tests with `test_package`

Run tests locally using e.g. `test_package('nimble', 'my-component')` to test `$REPO/packages/nimble/tests/testthat/test-my-component.R`.
Note that although can be using either `test_package()` or directly via `Rscript`, the `Rscript` method fails to detect some errors.

## Style

> Programs must be written for people to read,
> and only incidentally for machines to execute.
>
> --<cite>Harold Abelson, Structure and Interpretation of Computer Programs</cite>

Write code as if it is being published and will be read by your peers. (It is and will be.)
Ask yourself "Would I write this way in a preprint"?

Try to be consistent with the code you're modifying.
If in doubt, seek guidance from
[Google's R Style Guide](https://google.github.io/styleguide/Rguide.xml) or
[Hadley Wickham's style guide](http://adv-r.had.co.nz/Style.html).

For C++, install and use [`clang-format`](https://clang.llvm.org/docs/ClangFormat.html).
```sh
sudo apt-get install clang-format  # Linux.
brew install clang-format          # OS X.
```
Then after editing, let the robots tidy up
```sh
clang-format -i my_file.h my_file.cpp
```
or in the `$REPO/packages` directory, add your files to the `clang-format` target and simply
```sh
make clang-format
```

### Write clear comments

Write clear concise comments aimed at other developers.
If you write a great internal helper function,
and you want other developers to be able to use your code,
then write comment describing what your helper is and how to use it.
Write comments grammatically:
sentences should start with a capital letter and end with a period.
As yourself "Will my colleague be able to orient herself in this code with the aid of my comments"?

### Avoid submitting commented-out code

Avoid submitting commented-out code into the repo.
If you really want to remember the code, then remove it from the repo with a descriptive commit message.
If you really want the code to exist, then gate it by an `if()else` with an argument or global option.

To search git history for deleted code, use [`git grep`](https://git-scm.com/docs/git-grep). For example
```sh
$ git grep 'my_clever_sampler <- function' $(git log -g --pretty=format:%h)
```
You can even wrap this in the following shell script, say named `ggrep`:
```sh
#!/bin/sh
git grep "$@" $(git log -g --pretty=format:%h)
```

### Naming conventions

Try to be consistent.  We try to use `camelCase`.
```r
## Write grammatically correct properly punctuated comments.
## TODO(perrydv) This comment only makes sense for perrydv; ask him for details.

myVariable <- c(1, 2, 3)                               ## Variable names are lower case.
myFunction <- function(...) { return(list(...)) }      ## Function names are lower case.
MyRefClass <- setRefClass(Class = 'MyRefClass', ...)   ## Class names are Capitalized.
```
```c++
// Write grammatically correct properly punctuated comments.
// TODO(fritzo) This comment only makes sense to fritzo; ask him for details.

const double myVariable[] = {1, 2, 3};                 // Variable names are lower case.
void myFunction(int arg) { return arg; }               // Function names are lowere case.
class MyClass {                                        // Class names are upper case.
 public:
  void myMethod();
 private:                                              // Try to keep member variables private.
  myMemberVariable_;                                   // Member variables end with a trailing underscore _.
};
```

### Line width

In new code, try to keep line width to a maximum of **80 columns**.
This makes it much easier to use tooling:
- Side-by-side diffs are often easier to review than unified diffs, and prefer narrow code.
- `git blame` and `git grep` are clearer with narrow code.
- Software integration and debugging are easier when developers can view many files at once.
- Not everyone has an enormous screen: some developers respond to bugs using their phone.

## Debugging

Some Nimble errors are silenced by default.
While developing Nimble, you should enable more verbose error reporting.
To get verbose tracebacks, set the following (maybe in your `.Rprofile`)
```r
options(error = traceback)
options(show.error.locations = TRUE)
```
Rstudio is known to silence Nimble tracebacks.
The [recommmended workaround](https://support.rstudio.com/hc/en-us/community/posts/218745527/comments/223827227)
is to disable `Tools -> Global Options -> General -> [ ] Use debug error handler only when my code contains errors?`.

Make use of standard R debugging tools
```r
browser()                # Temporarily add this to code you are trying to debug.
debug(nimble:::foo)      # Call this before running something.
debugonce(nimble:::foo)  # Similar, but only debug on first call.
```
Nimble also has `nimbleOptions` that can be enabled for debugging:
- `debugRCfunProcessing`
- `debugNFProcessing`
- `debugSizeProcessing`
- `debugCppLineByLine`
