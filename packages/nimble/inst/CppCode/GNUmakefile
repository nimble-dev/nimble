SRC=RcppUtils.cpp Utils.cpp NamedObjects.cpp ModelClassUtils.cpp accessorClasses.cpp dists.cpp NFxnUtils.cpp Values.cpp nimDists.cpp nimbleGraph.cpp
OBJECTS=$(SRC:%.cpp=%.o)

ifndef R_SHARE_DIR
 R_SHARE_DIR=$(R_HOME)/share
endif

include $(R_HOME)/etc$(R_ARCH)/Makeconf

OS=$(shell uname)
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

%.o: %.cpp
	$(R_HOME)/bin/R CMD COMPILE $<


clean:
	-rm *.o libnimble$(DYLIB_EXT)
