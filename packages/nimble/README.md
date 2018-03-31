# UNIX

One can specify the path to where the Eigen header files are located
and also whether to link common C++ files into every DSO/DLL we create
or whether to treat the files as an extra library with:

```
R CMD INSTALL --configure-args='--with-eigen=/home/duncan/local --enable-lib' nimble
```
or, within R
```
install.packages("nimble", configure.args = c("--with-eigen=/home/duncan/local", "--enable-lib=true"), repos = NULL)
```


# Windows

Typically, you need the R developer tools (i.e., compiler, make, etc.) to use nimble.
Accordingly, it is quite straightforward to install the package from source as you will have the necessary tools
already installed. These are available from the [Rtools](https://cran.r-project.org/bin/windows/Rtools/) page on CRAN.

To install the package from source, from within R,
```r
install.packages("nimble", type = "source", INSTALL_opts = "--merge-multiarch")
```
or from a local copy of the source package,
```
install.packages("nimble_0.6-2.tar.gz", repos = NULL, INSTALL_opts = "--merge-multiarch")
```
Alternatively, use the shell command (in the DOS Command prompt)
```
R CMD INSTALL --merge-multiarch nimble_0.6-9.tar.gz
```
Of course, you can also compile directly from a clone of the git repository:
```
R CMD INSTALL --merge-multiarch nimble
```

The --merge-multiarch is necessary when using a version of R that supports both 32 and 64 bit.
This option to installation will ensure that  create both 32 and 64 bit installations.

## Creating a Windows Binary
```
R CMD build nimble
R CMD INSTALL --build --merge-multiarch nimble_0.6-9.tar.gz
```
We need to create the .tar.gz file first, hence the first command.
