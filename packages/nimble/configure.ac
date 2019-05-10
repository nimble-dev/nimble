dnl This is the configuration script for NIMBLE.
dnl Currently, it addresses having the header files from the Eigen software/library
dnl in non-standard locations on the machine.
dnl The idea is to generate a Makevars file in inst/make that will be copied or linked
dnl to the user's current directory or the directory in which the C++ code NIMBLE creates.


AC_INIT(DESCRIPTION)

dnl We should probably get the C++ compiler from R's configuration script to ensure compatability.

NEED_MAKEVARS_FILE=TRUE

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
    AC_MSG_ERROR([Could not determine R_HOME.])   
fi

dnl # Define -fPIC rather than -fpic when possibly, to overcome -fpic's limit on
dnl # the number of exported symbols.
dnl CXXPICFLAGS="-fPIC"
dnl case "${CXX}" in
dnl   ## Intel compilers
dnl   *icpc|*icc)
dnl     CXXPICFLAGS="-fpic"
dnl     ;;
dnl   ## Portland Group
dnl   *pgCC)
dnl     CXXPICFLAGS="-fpic"
dnl     ;;
dnl esac

dnl The caller can specify where their Eigen/ directory is using, e.g.,
dnl    --configure-args='--with-eigen=/home/bob/local'
dnl or by setting the environment variable EIGEN_DIR
AC_ARG_WITH(eigen, [  --with-eigen            path for the directory containing Eigen/ directory with all of the #include files], EIGEN_DIR="${withval}")

dnl Similarly, caller can specify their CppAD directory via
dnl    --configure-args='--with-cppad=/home/bob/local/cppad'
dnl or by setting the environment variable CPPAD_DIR
AC_ARG_WITH(cppad, [  --with-cppad            path for the directory containing cppad/ directory], CPPAD_DIR="${withval}")

AC_ARG_ENABLE(dylib, [  --enable-dylib            use a single shared libnimble rather than linking code into each DLL/DSO generated when compiling generated code],
                         ENABLE_LIB=${enableval},
if ! test `uname` = "Linux" ; then
   ENABLE_LIB=false
else
   ENABLE_LIB=true
fi
)

if test -n "${NO_LIBNIMBLE}" ; then
  ENABLE_LIB=false
  echo "Disabling libnimble."
fi

if test -n "${LIBNIMBLE}" ; then
  ENABLE_LIB=true
  echo "Enabling libnimble."
fi

dnl --enable-debugging or --enable-debugging=yes will keep any -g flags from R's Makeconf
dnl default (omission of --enable-debugging) will invoke use of customizeMakeconf in inst/CppCode
dnl    which in turn makes a local copy of Makeconf and replaces "-g " with " ".
AC_ARG_ENABLE(debugging, [  --enable-debugging            compiled libnimble with debugging information],
ENABLE_DEBUGGING=${enableval},
ENABLE_DEBUGGING=no)

if test "${ENABLE_DEBUGGING}" = "yes"; then
   echo "enabling R's default debugging information in libnimble"
   GNUMAKEFILE_INCLUDE_MAKECONF="\$(R_HOME)/etc/\$(R_ARCH)/Makeconf"
else
   GNUMAKEFILE_INCLUDE_MAKECONF="Makeconf"
fi
AC_SUBST(GNUMAKEFILE_INCLUDE_MAKECONF)

if test -n "$EIGEN_DIR" ; then
  echo "Checking eigen directory $EIGEN_DIR"
  if ! test -d "$EIGEN_DIR" ; then
     echo "$EIGEN_DIR is not a directory"
  fi

  if ! test -d "$EIGEN_DIR/Eigen" ; then
     echo "$EIGEN_DIR does not contain the directory Eigen"
     exit 1;
  fi

  if ! test -f "$EIGEN_DIR/Eigen/Dense" ; then
     echo "$EIGEN_DIR/Eigen does not contain the file Dense"
     exit 1;
  fi

  if ! test "${EIGEN_DIR}" = "" ; then
    EIGEN_INC="-I${EIGEN_DIR}"
  fi
fi

if test -n "$CPPAD_DIR" ; then
  echo "Checking cppad directory CPPAD_DIR"
  if ! test -d "$CPPAD_DIR" ; then
     echo "$CPPAD_DIR is not a directory"
  fi

  if ! test -d "$CPPAD_DIR/cppad" ; then
     echo "$CPPAD_DIR does not contain the directory cppad"
     exit 1;
  fi

  if ! test -f "$CPPAD_DIR/cppad/configure.hpp" ; then
     echo "$CPPAD_DIR/cppad does not contain the file configure.hpp"
     exit 1;
  fi
fi

if test -n "CPPAD_DIR" && ! test "${CPPAD_DIR}" = "" ; then
   echo "Setting CPPFLAGS to find CppAD"
   export CPPFLAGS="-I${CPPAD_DIR}"
fi

if test -n "$EIGEN_DIR" && ! test "${EIGEN_DIR}" = "" ; then
   echo "Setting CPPFLAGS to find Eigen"
   export CPPFLAGS="-I${EIGEN_DIR}"
fi


dnl Next 5 lines are from from https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Using-C_002b_002b11-code, but with CXX11 modified to CXX1X per output of R CMD config
dnl 2018-07-25: now modified to CXX11 as that seems to be standard on more recent R versions
CXX11=`"${R_HOME}/bin/R" CMD config CXX11`
CXX11STD=`"${R_HOME}/bin/R" CMD config CXX11STD`
CXX="${CXX11} ${CXX11STD}"
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXX11FLAGS`
AC_LANG(C++)
dnl 2018-07-25: moved next line from early in this file because it seems to trigger various configure checks based on g++ even when R is configured to use clang; see issue 766
dnl however note that having it here seems extraneous as it produces same configure file without it
AC_PROG_CXX  

dnl Next 3 lines are from old config.ac, prior to turning on C++11 and using above 5 lines
dnl CXX=`"${R_HOME}/bin/R" CMD config CXX`
dnl CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
dnl AC_LANG([C++])

echo "working on eigen"

AC_CHECK_HEADERS([Eigen/Dense],
                 [echo "Found Eigen/Dense";
dnl                rm -rf inst/Eigen_local
		  if test -e "NOT_IN_TARBALL" ; then
		   echo "Providing your own Eigen only works if installing from a tarball.  Try R CMD BUILD nimble; R CMD INSTALL --configure-args='--with-eigen=YOUR_EIGEN_PATH' nimble_0.6-1.tar.gz or something similar."
		  else
		   echo "Not installing Eigen with NIMBLE"
		   rm -r inst/include/Eigen
		  fi
                 ],
                 [echo "Using the version of Eigen provided with NIMBLE";
		  EIGEN_DIR=""
                  ])

dnl Somewhat confused if this is redundant with the CPPFLAGS entry above
if test -n "$CPPAD_DIR" ; then
echo "Using CppAD directory ${CPPAD_DIR}"
CPPAD_INC="-I$CPPAD_DIR"
fi

if test -n "$EIGEN_DIR" ; then
echo "Using Eigen directory ${EIGEN_DIR}"
EIGEN_INC="-I$EIGEN_DIR"
fi

dnl AC_SUBST(CXXPICFLAGS)
AC_SUBST(CPPAD_INC)
AC_SUBST(EIGEN_INC)
AC_SUBST(NEED_MAKEVARS_FILE)
AC_SUBST(R_PACKAGE_DIR)

dnl As of staged installation (2019-03-01)
dnl NIMBLE_DIR, NIMBLE_INC_DIR, NIMBLE_LIB_DIR, RPATH for
dnl default Linux enable-lib=true configuration are only
dnl being set temporarily to template values (__FOO__)
dnl at installation and then set at run-time via createMakevars().

if test "${ENABLE_LIB}" != "false" ; then
 NIMBLE_DIR="__NIMBLE_DIR__"
else
 NIMBLE_DIR="WILL_BE_SET_IN_LOCAL_MAKEVARS"
fi
AC_SUBST(NIMBLE_DIR)

if test "${ENABLE_LIB}" != "false" ; then
 NIMBLE_INC_DIR="__NIMBLE_INC_DIR__"
else
 NIMBLE_INC_DIR="WILL_BE_SET_IN_LOCAL_MAKEVARS"
fi
AC_SUBST(NIMBLE_INC_DIR)

if test "${ENABLE_LIB}" = "yes" ; then
 ENABLE_LIB="TRUE"
elif test "${ENABLE_LIB}" = "no" ; then
  ENABLE_LIB="FALSE"
fi

if test "${ENABLE_LIB}" != "false" ; then
 NIMBLE_LIB_DIR="__NIMBLE_LIB_DIR__"
else
 NIMBLE_LIB_DIR="WILL_BE_SET_IN_LOCAL_MAKEVARS"
fi
AC_SUBST(NIMBLE_LIB_DIR)


R_ENABLE_LIB=`echo "${ENABLE_LIB}" | tr "[[a-z]]" "[[A-Z]]"`
echo "ENABLE_LIB=${ENABLE_LIB}"
AC_SUBST(R_ENABLE_LIB)

dnl Note RPATH only now used at run-time via make/Makevars*
if ! test `uname` = "Darwin" ; then
RPATH="__RPATH__"   ## -Wl,-rpath=${NIMBLE_DIR}/CppCode"
fi
AC_SUBST(RPATH)

dnl if test "${ENABLE_LIB}" != "false" ; then
dnl  #XX
dnl  PKG_LIBS="-L../inst/CppCode -lnimble $RPATH"
dnl #  PKG_LIBS="-L. -lnimble $RPATH"
dnl fi
dnl AC_SUBST(PKG_LIBS)

USE_REGISTRATION=".NimbleUseRegistration = TRUE"
AC_SUBST(USE_REGISTRATION)

AC_OUTPUT(inst/make/Makevars inst/make/Makevars_lib R/config.R src/Makevars)

if test "${ENABLE_DEBUGGING}" != "yes"; then
  echo "removing debugging flags"
  (cd inst/CppCode ; ${MAKE} -f GNUmakefile_customizeMakeconf localMakeconf)
fi

cp inst/CppCode/dists.cpp src
cp inst/CppCode/RcppUtils.cpp src
cp inst/CppCode/nimbleGraph.cpp src
cp inst/CppCode/eigenUsingClasses.cpp src

if test "${ENABLE_LIB}" != "false" ; then
  echo "Building libnimble"
  (cd inst/CppCode ; ${MAKE} -f${GNUMAKEFILE_INCLUDE_MAKECONF} -fGNUmakefile)

else
  echo "Building libnimble.a"
  (cd inst/CppCode ; ${MAKE} -f${GNUMAKEFILE_INCLUDE_MAKECONF} -fGNUmakefile libnimble.a)
fi

if test "${ENABLE_DEBUGGING}" != "yes"; then
  (cd inst/CppCode ; ${MAKE} -f GNUmakefile_customizeMakeconf clean)
fi
