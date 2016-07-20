Classic BUGS examples
---------------------

Some of the BUGS examples have been modified to run with JAGS, and have
been turned into a test suite.

Each subdirectory contains a number of JAGS script files named test1.cmd,
test2.cmd, ... etc. designed to run the examples as they were presented
in the classic BUGS examples book.

There are also R script files test1.R, test2.R ... to run the
examples using the rjags interface from JAGS to R.

Some of the examples have been modified. In such cases, there is a
ReadMe file explaining the changes. Some examples cannot be run at all,
and these are not included in the testing process.

Make targets
------------

The test suite uses GNU make.  To run the tests, change directory into
either "vol1" or "vol2" and invoke make with one of the following targets.

make check  Will run all the examples and check the results against
            the benchmark run.  The check ensures that the posterior
            mean of the monitored nodes is not more than 0.15 standard
            deviations away from the benchmark mean. 

make Rcheck Runs examples using the rjags interface.

make clean  Cleans up any files that may have been created during
            make check, or during debugging

make bench  Creates new benchmark results

make distclean Deletes the benchmark results

In order to run make for a specific example, use the variable "EXAMPLES":
make check EXAMPLES=line

Debugging an example
--------------------

If "Make check" fails for example foo, then make will report an error in
make taget "check_foo".  The log file "check.log" in subdirectory "foo"
can be inspected to see what went wrong. This log file includes the
output produced by JAGS when running the model and by R when checking
the output. Similarly, "make Rcheck" exands to targets of the form
"Rcheck_foo" and the log file is "Rcheck.log".

Note that the tests in "make check" and "make Rcheck" are not deterministic and 
will occasionally fail just by chance.  You may run "make check" again 
when this happens to see if the error is consistently reproduced.

Re-running the examples
-----------------------

Once all the tests have been successfully completed in a sub-directory,
an empty marker file will be created (check.OK or Rcheck.OK) and the 
corresponding checks will be skipped.  To re-run the checks, first use
"make clean" to remove the marker files.
 
Parallel make
-------------

The examples in separate directories can be run in parallel for increased
efficiency.  Use the "-j" argument to run examples in parallel, e.g

make -j4 check

Will run 4 parallel make jobs.
