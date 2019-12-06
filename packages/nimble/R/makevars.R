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

  ## This old condition evaluates to TRUE even when we should be generating a new makevars
  ## For now I am simplifying by conditioning only on .useLib
  ##  haveContents = !missing(pkgFlags) || !missing(pkgLibs) || length(list(...))

  ##  if(!haveContents) {
  if(.useLib) {
     if(!file.exists(.copyFrom))
         stop("No default Makevars file")

     ## file.copy(.copyFrom, target)  # file.link won't work across file systems.
     contents <- readLines(.copyFrom)
     NIMBLE_INC_DIR <- normalizePath(system.file("include", package = "nimble"), winslash = "\\", mustWork=FALSE)
     NIMBLE_LIB_DIR <- normalizePath(system.file("CppCode", package = "nimble"), winslash = "\\", mustWork=FALSE)
     NIMBLE_DIR <- normalizePath(system.file(package = "nimble"), winslash = "\\", mustWork=FALSE)
     RPATH <- sprintf("-Wl,-rpath %s", normalizePath(system.file("CppCode", package = "nimble"), winslash = "\\", mustWork=FALSE))
     contents <- gsub("__NIMBLE_INC_DIR__", NIMBLE_INC_DIR, contents)
     ## NIMBLE_LIB_DIR is being set relative to NIMBLE_DIR so this is not really necessary,
     ## but keeping in case of unexpected corner case.
     contents <- gsub("__NIMBLE_LIB_DIR__", NIMBLE_LIB_DIR, contents)
     contents <- gsub("__NIMBLE_DIR__", NIMBLE_DIR, contents)
     contents <- gsub("__RPATH__", RPATH, contents)
     cat(contents, file = target, sep = "\n")
     return(target)
  }

  args = list(...)
  if(!missing(pkgFlags))
     args$PKG_CPPFLAGS = pkgFlags
  if(!missing(pkgLibs))
     args$PKG_LIBS = pkgLibs
  args = sapply(args, as.character)
#  cat(sprintf("%s=%s", names(args), args), sep = "\n", file = target)
  genLocalMakevars(target, args, .useLib)
  
  target
}

genLocalMakevars =
function(target, vars = character(), .useLib = UseLibraryMakevars)
{
##    cat("creating local makeVars in", target, "\n")
    inc.make = system.file("make", if(.useLib) 
                                     "Makevars_lib" 
                                   else if(.Platform$OS.type == "windows")
                                     "Makevars.win"
                                   else
                                     "Makevars", package = "nimble")

    cppad_inc <- paste0("'-I\"", nimbleOptions('CppADdir'), "\"'")
    vars = c(EIGEN_INC = "", ## AutoconfInfo$eigenInc, ## we used to generate an AutoconfInfo list. We'll need a new mechanism if a local Makevars needs to be generated and the user has non-nimble-provided Eigen
             CPPAD_INC = "", ##cppad_inc,  ## Currently, require cppad folder to be placed within include folder.  Todo: get CppAD_directory functional.
             NIMBLE_INC_DIR =  normalizePath(system.file("include", package = "nimble"), winslash = "\\", mustWork=FALSE),
             NIMBLE_LIB_DIR =  normalizePath(system.file("CppCode", package = "nimble"), winslash = "\\", mustWork=FALSE),
             NIMBLE_DIR =  normalizePath(system.file(package = "nimble"), winslash = "\\", mustWork=FALSE),
             RPATH = sprintf("-rpath %s", normalizePath(system.file("CppCode", package = "nimble"), winslash = "\\", mustWork=FALSE)),
             vars)
    varDefs = mapply(function(id, val) paste(id, val, sep = "="), names(vars), vars)

       # replace any spaces in the path with \<space>, but need \\\\ to get the single \ 
    inc.make = gsub(" ", "\\\\ ", inc.make)
    content = c(varDefs, "", sprintf('include %s', inc.make))

    cat(content, file = target, sep = "\n")
}




sameDir =
function(a, b)
{
    path.expand(a) == path.expand(b)
}
