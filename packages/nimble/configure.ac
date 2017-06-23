dnl This is the configuration script for NIMBLE.
dnl Currently, it addresses having the header files from the Eigen software/library
dnl in non-standard locations on the machine.
dnl The idea is to generate a Makevars file in inst/make that will be copied or linked
dnl to the user's current directory or the directory in which the C++ code NIMBLE creates.


AC_INIT(DESCRIPTION)

dnl We should probably get the C++ compiler from R's configuration script to ensure compatability.
AC_PROG_CXX

NEED_MAKEVARS_FILE=TRUE

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


CXX=`"${R_HOME}/bin/R" CMD config CXX`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
AC_LANG([C++])
dnl  AC_LANG_PUSH([C++]);

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
if test -n "CPPAD_DIR" ; then
echo "Using CppAD directory ${CPPAD_DIR}"
CPPAD_INC="-I$CPPAD_DIR"
fi

if test -n "$EIGEN_DIR" ; then
echo "Using Eigen directory ${EIGEN_DIR}"
EIGEN_INC="-I$EIGEN_DIR"
fi

AC_SUBST(CPPAD_INC)
AC_SUBST(EIGEN_INC)
AC_SUBST(NEED_MAKEVARS_FILE)
AC_SUBST(R_PACKAGE_DIR)

if test "${ENABLE_LIB}" != "false" ; then
 NIMBLE_DIR=${R_PACKAGE_DIR}
 echo "NIMBLE_DIR = ${NIMBLE_DIR}"
else
 NIMBLE_DIR="WILL_BE_SET_IN_LOCAL_MAKEVARS"
fi
AC_SUBST(NIMBLE_DIR)

if test "${ENABLE_LIB}" != "false" ; then
 NIMBLE_INC_DIR=${NIMBLE_DIR}/include
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
 NIMBLE_LIB_DIR=${NIMBLE_DIR}/CppCode
else
 NIMBLE_LIB_DIR="WILL_BE_SET_IN_LOCAL_MAKEVARS"
fi
AC_SUBST(NIMBLE_LIB_DIR)


R_ENABLE_LIB=`echo "${ENABLE_LIB}" | tr "[[a-z]]" "[[A-Z]]"`
echo "ENABLE_LIB=${ENABLE_LIB}"
AC_SUBST(R_ENABLE_LIB)

if ! test `uname` = "Darwin" ; then
RPATH="-Wl,-rpath=${NIMBLE_DIR}/CppCode"
fi
AC_SUBST(RPATH)

if test "${ENABLE_LIB}" != "false" ; then
  #XX
  PKG_LIBS="-L../inst/CppCode -lnimble $RPATH"
#  PKG_LIBS="-L. -lnimble $RPATH"
fi
AC_SUBST(PKG_LIBS)

USE_REGISTRATION=".NimbleUseRegistration = TRUE"
AC_SUBST(USE_REGISTRATION)

AC_OUTPUT(inst/make/Makevars inst/make/Makevars_lib R/config.R src/Makevars)



if test "${ENABLE_LIB}" != "false" ; then
 echo "Building libnimble"
 (cd inst/CppCode ; ${MAKE} )
# cp inst/CppCode/libnimble.dylib src
# if ! test -e $R_PACKAGE_DIR/libs/ ; then
#    mkdir -p $R_PACKAGE_DIR/libs/
# fi
 #cp inst/CppCode/libnimble.dylib $R_PACKAGE_DIR/libs/libnimble.dylib 
else 
#  cp inst/CppCode/*.cpp src
  cp inst/CppCode/dists.cpp src
  cp inst/CppCode/RcppUtils.cpp src
  cp inst/CppCode/nimbleGraph.cpp src
  cp inst/CppCode/eigenUsingClasses.cpp src

  (cd inst/CppCode ; ${MAKE} libnimble.a)
fi
