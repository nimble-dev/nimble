.onAttach <- function(libname, pkgname) {
    release <- utils::packageDescription("nimble", field = "Version")
    #date <- utils::packageDescription("nimble", field = "Date")
    packageStartupMessage("nimble version ", release, " is loaded.",
                          "\nFor more information on NIMBLE and a User Manual,",
                          "\nplease visit https://R-nimble.org.",
                          "\n\nNote for advanced users who have written their own MCMC samplers:",
                          "\n  As of version 0.13.0, NIMBLE's protocol for handling posterior",
                          "\n  predictive nodes has changed in a way that could affect user-defined",
                          "\n  samplers in some situations. Please see Section 15.5.1 of the User Manual."
                          )
}

if(!exists("NeedMakevarsFile"))
 NeedMakevarsFile = FALSE
if(!exists("UseLibraryMakevars"))
  UseLibraryMakevars = FALSE
