createMakevars =
#
# Create a Makevars file in the specified dir(ectory).
# This either copies a Makevars file from the nimble installation
# or creates one with the contents specified via the pkgFlags, pkgLibs, ....
#
# There are two Makevars files in nimble.  One compiles the .cpp files in inst/include (currently)
# used in all DSOs we create and links them with the generated code to create a single .so.  The other combines this 
# common code into their own library and then links against this. There are issues with being able to find
# that DSO/library at run time. We make some efforts to put the path into the linking step.
# The goal of this second approach is to allow the code to be shared across all DSOs we load
# rather than repeat it each time. So this is an effort to save space.
#
function(pkgFlags, pkgLibs, ..., dir = getwd(), 
         .copyFrom = getOption("nimble.Makevars.file",
                               system.file("make", 
                                           sprintf("%s%s",
                                                   if(.useLib) 
                                                      "Makevars_lib" 
                                                   else 
                                                      "Makevars", 
                                                   if(.Platform$OS.type == "windows") ".win" else ""),
                                                   package = "nimble")), 
         .force = FALSE, .useLib = UseLibraryMakevars)
{
  target = sprintf("%s%s%s%s", dir, .Platform$file.sep, "Makevars", if(.Platform$OS.type == "windows") ".win" else "")

  if(file.exists(target) && !.force)
     stop(paste(target, "already exists"))

  haveContents = !missing(pkgFlags) || !missing(pkgLibs) || length(list(...))

  if(!haveContents) {
     if(!file.exists(.copyFrom))
         stop("No default Makevars file")

     file.copy(.copyFrom, target)  # file.link won't work across file systems.
     return(target)
  }

  args = list(...)
  if(!missing(pkgFlags))
     args$PKG_CPPFLAGS = pkgFlags
  if(!missing(pkgLibs))
     args$PKG_LIBS = pkgLibs
  args = sapply(args, as.character)
  cat(sprintf("%s=%s", names(args), args), sep = "\n", file = target)

  target
}
