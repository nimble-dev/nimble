.onLoad <- function(libname, pkgname) {
    packageStartupMessage('For more information on NIMBLE and a User Manual, please visit http://R-nimble.org')
}

if(!exists("NeedMakevarsFile"))
 NeedMakevarsFile = FALSE
if(!exists("UseLibraryMakevars"))
  UseLibraryMakevars = FALSE
