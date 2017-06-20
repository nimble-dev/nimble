build: nimble/configure
	@echo "Building nimble package in this directory"
	(R CMD build nimble)

configure: nimble/configure.ac
	(cd nimble; autoconf)

clang-format:
	clang-format -style=file \
	  nimble/inst/include/nimble/NimArrBase.h \
	  nimble/inst/include/nimble/NimArr.h \
	  nimble/inst/include/nimble/predefinedNimbleLists.h \
	  nimble/inst/CppCode/predefinedNimbleLists.cpp \
	  -i
