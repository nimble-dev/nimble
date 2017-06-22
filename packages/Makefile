.PHONY: FORCE

build: nimble/configure
	@echo "Building nimble package in-place"
	(R CMD build nimble)

configure: nimble/configure.ac
	(cd nimble; autoconf)

clang-format: FORCE
	clang-format -style=file \
	  nimble/inst/include/nimble/NimArrBase.h \
	  nimble/inst/include/nimble/NimArr.h \
	  nimble/inst/include/nimble/predefinedNimbleLists.h \
	  nimble/inst/CppCode/predefinedNimbleLists.cpp \
	  -i

test: FORCE
	cd .. ; ./run_tests.R --parallel

clean: FORCE
	git clean -dfX

FORCE:

