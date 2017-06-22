.PHONY: FORCE

VERSION=$(shell grep Version: nimble/DESCRIPTION | grep -o '[^ ]*$$')

version: FORCE
	@echo $(VERSION)

configure: nimble/configure.ac
	(cd nimble; autoconf)

nimble_$(VERSION).tar.gz: nimble/configure FORCE
	@echo "Building nimble package in this directory"
	(R CMD build nimble)

build: nimble_$(VERSION).tar.gz

check: build FORCE
	R CMD check --as-cran nimble_$(VERSION).tar.gz

test: FORCE
	cd .. ; ./run_tests.R --parallel

clang-format: FORCE
	clang-format -style=file \
	  nimble/inst/include/nimble/NimArrBase.h \
	  nimble/inst/include/nimble/NimArr.h \
	  nimble/inst/include/nimble/predefinedNimbleLists.h \
	  nimble/inst/CppCode/predefinedNimbleLists.cpp \
	  -i

clean: FORCE
	git clean -dfX

FORCE:

