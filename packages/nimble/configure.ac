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


if test -n "$EIGEN_DIR" && ! test "${EIGEN_DIR}" = "" ; then
   echo "Setting CPPFLAGS to find Eigen"
   export CPPFLAGS="-I${EIGEN_DIR}"
fi

AC_LANG([C++])
dnl  AC_LANG_PUSH([C++]);

AC_CHECK_HEADERS([Eigen/Dense], 
                 [echo "Found Eigen/Dense";
dnl                rm -rf inst/Eigen_local
                 ], 
                 [echo "Cannot find Eigen/Dense include file.";
                  echo "Using the version of Eigen provided locally with the R package";
                  echo "To avoid this, re-install and specify the path containing the Eigen/ directory via";
                  echo "   --with-eigen (see INSTALL for syntax) or set the environment variable EIGEN_DIR";
		  if test -e "NOT_IN_TARBALL" ; then
		   echo "found NOT_IN_TARBALL"
		   echo `pwd`
                   cp -r inst/Eigen_local inst/include/Eigen 
		  else
		   echo "did not find NOT_IN_TARBALL"
		   echo `pwd`
		   mv inst/Eigen_local inst/include/Eigen 
		  fi
                  EIGEN_DIR=""
                  ])


if test -n "$EIGEN_DIR" ; then
echo "Using Eigen directory ${EIGEN_DIR}"
EIGEN_INC="-I$EIGEN_DIR"
fi


AC_SUBST(EIGEN_INC)
AC_SUBST(NEED_MAKEVARS_FILE)
AC_SUBST(R_PACKAGE_DIR)

NIMBLE_DIR=${R_PACKAGE_DIR}
echo "NIMBLE_DIR = ${NIMBLE_DIR}"
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


R_ENABLE_LIB=`echo "${ENABLE_LIB}" | tr "[a-z]" "[A-Z]"`
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
  cp inst/CppCode/*.cpp src
  (cd inst/CppCode ; ${MAKE} libnimble.a)
fi
