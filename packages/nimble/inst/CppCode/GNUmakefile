SRC=RcppUtils.cpp EigenTypedefs.cpp eigenUsingClasses.cpp RcppNimbleUtils.cpp Utils.cpp NamedObjects.cpp ModelClassUtils.cpp accessorClasses.cpp dists.cpp nimDists.cpp dllFinalizer.cpp nimbleGraph.cpp smartPtrs.cpp predefinedNimbleLists.cpp nimOptim.cpp
OBJECTS=$(SRC:%.cpp=%.o)

# nimbleGraph.cpp is only needed in SRC for --enable-dylib=true, but it's not a big deal to leave include it either way

ifndef R_SHARE_DIR
 R_SHARE_DIR=$(R_HOME)/share
endif

# The following two lines attempt to reduce the size of libnimble.a or equivalent by removing -g flags.
# It is OK if the R script silently fails; in that case Makeconf will simply not be found.
# Note that we don't precede this with a fallback `-include $(R_HOME)/etc$(R_ARCH)/Makeconf`,
# since that causes "make variable override" warnings on OS X.
##$(shell $(R_HOME)/bin$(R_ARCH)/Rscript --vanilla ./customizeMakeconf.R)
##-include ./Makeconf

## reverting back to this include of Makeconf
-include $(R_HOME)/etc$(R_ARCH)/Makeconf
include ./Makevars

ifndef OS
 OS=$(shell uname)
endif

ifeq ($(OS), Darwin)
INSTALL_NAME=-install_name "$(R_PACKAGE_DIR)/CppCode/libnimble$(DYLIB_EXT)"
#Relative paths don't work on OSX. "../CppCode/libnimble$(DYLIB_EXT)"
endif

libnimble$(DYLIB_EXT): $(OBJECTS)
	@echo "linking libnimble$(DYLIB_EXT)"
	@echo "$(OS)   $(INSTALL_NAME)"
	$(CXX) $(SHLIB_CXXLDFLAGS) -o $@  $(OBJECTS) $(INSTALL_NAME)
	-rm *.o

#	cp libnimble$(DYLIB_EXT) ../../src

#%.o: %.cpp
#	$(R_HOME)/bin$(R_ARCH)/Rcmd COMPILE $<
# Was R CMD COMPILE

clean:
	-rm -f $(OBJECTS) libnimble.$(DYLIB_EXT) libnimble.dylib libnimble$(ARCH_SUFFIX).a
# DYLIB_EXT does not seem to be defined unless this is called within R CMD.  e.g. from ./cleanup, DYLIB_EXT is not defined

# libnimble32.a libnimble64.a
libnimble.a  libnimble_x64.a libnimble_i386.a: $(OBJECTS)
	echo "Making nimble $@"      
	echo "$(OBJECTS)"
	$(AR) rcs $@ $(OBJECTS)
	-rm $(OBJECTS)

