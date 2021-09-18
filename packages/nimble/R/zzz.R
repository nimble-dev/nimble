.onAttach <- function(libname, pkgname) {
    release <- utils::packageDescription("nimble", field = "Version")
    #date <- utils::packageDescription("nimble", field = "Date")
    packageStartupMessage("nimble version ", release, " is loaded.",
                          "\nFor more information on NIMBLE and a User Manual,",
                          "\nplease visit https://R-nimble.org.")
}

if(!exists("NeedMakevarsFile"))
 NeedMakevarsFile = FALSE
if(!exists("UseLibraryMakevars"))
  UseLibraryMakevars = FALSE
