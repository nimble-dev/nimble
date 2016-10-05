login as paciorek with s1x2 as password
scp paciorek@smeagol.berkeley.edu:~/research/perry/solaris/pkg.oracle.com.key.pem /tmp/.
scp paciorek@smeagol.berkeley.edu:~/research/perry/solaris/pkg.oracle.com.certificate.pem /tmp/.
su - # s1x2
pkg set-publisher \
       	    -k /tmp/pkg.oracle.com.key.pem \
       	    -c /tmp/pkg.oracle.com.certificate.pem \
       	    -G "*" -g https://pkg.oracle.com/solarisstudio/release solarisstudio

pkg publisher solarisstudio | grep Mirror

If the output is empty you are all set. If not remove unrelated mirrors by running:
	$ sudo pkg set-publisher -M http://mirror1.x.com -M http://mirror2.y.com ... solarisstudio
	                
pkg install solarisstudio-123
# puts in /opt

export PATH=/opt/solarisstudio12.3/bin:$PATH

solaris 'add more software'; search for iconv and select iconv/extra, unicode, unicode-core, utf-8, text/locale
# didn't work - R manual says need GNU

wget http://ftp.gnu.org/pub/gnu/libiconv/libiconv-1.14.tar.gz
./configure
make 
make install

install solaris xz package (not new enough)

# this didn't work
pkgadd -d http://get.opencsw.org/now
/opt/csw/bin/pkgutil -U
/opt/csw/bin/pkgutil -y -i xz 
/usr/sbin/pkgchk -L CSWxz # list files
# hmm that install of xz may have done iconv too
export PATH=/opt/csw/bin:$PATH
export LD_LIBRARY_PATH=/opt/csw/lib

# this works
wget http://tukaani.org/xz/xz-5.2.2.tar.gz
configure/make/makeinstall

# issue with finding java binary (perhaps issue with stringi, hoepfully)
export PATH=/usr/java/bin:$PATH

/opt/csw/bin/pkgutil -a g++
/opt/csw/bin/pkgutil -y -i gcc4g++
/opt/csw/bin/pkgutil -y -i gcc4gfortran

mkdir /usr/lib/R
wget https://cran.r-project.org/src/base/R-3/R-3.3.1.tar.gz
tar -xvzf R-3.3.1.tar.gz
# put config.site from CRAN into R directory (otherwise get -KPIC in g++ calls, and that fails)
CC="/opt/csw//bin/gcc"
CFLAGS=-O2
CPPFLAGS="-I/opt/csw/include -I/usr/local/include"
F77="/opt/csw//bin/gfortran"
FFLAGS=-O2
CXX="/opt/csw//bin/g++"
CXXFLAGS=-O2
FC=$F77
FCFLAGS=-O2
LDFLAGS="-L/opt/csw/lib -L/usr/local/lib"

cd R-3.3.1
./configure MAKE=gmake  # try this to avoid having to hack gmake below
make

mv /usr/bin/awk{,-}
ln -s /usr/bin/nawk /usr/bin/awk # some issue with awk and igraph
install.packages(c('igraph','coda'), repos = 'https://cran.cnr.berkeley.edu')

# use this to fix error with igraph install
https://github.com/igraph/rigraph/pull/128
# in particular this patch to src/uuid/gen_uuid.c
-#ifdef SIOCGIFHWADDR
+#if defined(SIOCGIFHWADDR) && (!defined(__sun__))
https://github.com/igraph/rigraph/pull/128/files

# didn't bother with make install
bin/R
install.packages("nimble", repos = 'https://cran.cnr.berkeley.edu')

# issue with stringi so remove dependnecies

# then I get an issue with make
# so symlink /usr/bin/make to /usr/bin/gmake since CRAN uses gmake somehow

export _R_CHECK_FORCE_SUGGESTS_=0  # so doesn't need various Suggests including Linux JAGS
/usr/lib/R/R-3.3.1/bin/R CMD check --as-cran nimble_0.6-1.tar.gz
