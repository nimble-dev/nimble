
One can specify the path to where the Eigen header files are located
and also whether to link common C++ files into every DSO/DLL we create
or whether to treat the files as an extra library with:

R CMD INSTALL . --configure-args='--with-eigen=/home/duncan/local --enable-lib'

install.packages(".", configure.args = c("--with-eigen=/home/duncan/local", "--enable-lib=true"), repos = NULL)
