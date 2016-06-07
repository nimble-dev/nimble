#
# The code in this file is to deal with a binary installation.
# <aside>I'll say it again. This is a lot of work for no good reason. People installing
#  nimble who want to use it are overwhelmingly going to use it to compile code generated
#  by nimble and so they will need to have the developer tools installed. So they could just install
#  the package from source and none of this would be a problem.
#  Similarly, we can fix this with a simple change to R, i.e., a script that gets run when we install a binary
#  package and then we can change settings (compile, link and run-time locations). (I should have committed this a long time ago - mea culpa).
#</aside>


.onLoad =
    # We don't have to do this on load, but just at the first time the library is needed.
function(...)
{
   if(!sameDir(system.file(package = "nimble"), AutoconfInfo$rpkg_install_dir)) {
       message("Package was moved after installation")
       NimbleSessionInfo$soDir = system.file("CppCode", package = "nimble")
       NimbleSessionInfo$pkgMoved = TRUE
       if(FALSE && UseLibraryMakevars && Sys.info()["sysname"] == "Darwin")
               # See if we have the write permissions to change the id in place.
          tryCatch(rebakeSOPath(libnimbleDLL()),
                     error = function(e) {
                             # If we don't have permission to rewrite the dylib, then
                            copyNimbleDylib()
                         })
   } else
       NimbleSessionInfo$pkgMoved = FALSE
}


copyNimbleDylib =
    #
    # create a copy of the libnimble.dylib file in another directory.
    # The expectation is that is called for a binary installation on OSX
    # and where we don't have permission to write to the libnimble.dylib file
    # e.g., this is a centrally installed package shared by multiple users.
    #
    #
function(dir = getOption("LibNimbleDir", tempdir()))
{
    NimbleSessionInfo$soDir = dir
    f = sprintf("%s/%s", dir, basename(libnimbleDLL()))
    cat("Copying .so\n")       
    file.copy(libnimbleDLL(), dir)
    Sys.chmod(f, "0755")
    # Change the permissions on the file.
    cat("Setting new id in .dylib\n")
    rebakeSOPath(f)
    f
}


libnimbleDLL =
function()
{
  list.files(system.file("CppCode", package = "nimble"), pattern = "(so|dll|dylib)$", full.names = TRUE)
}

rebakeSOPath =
function(filename)
{

   # Should check to see if we need to change the id.
  cat("rebaking")
  # /Users/duncan/RTestPackages/nimble/CppCode/libnimble.dylib  /Users/duncan/RTestPackages/nimble/CppCode/libnimble.dylib
  # XXX needs install_name_tool to be installed on the machine (not necessarily the machine on which the *package* source was compiled)
  # This is part of CommandLineTools
  # /Library/Developer/CommandLineTools/usr/bin/install_name_tool
  # We may need to locate it explicitly.

  #
  # Should check the length of the name
  #

# print( getCurrentSOName(filename) )
#  cmd = sprintf("install_name_tool -id  %s %s", filename, filename)
#  system(cmd, intern = TRUE)  
#  err = textConnection("e", open = "w", local = TRUE)
#  out = textConnection("o", open = "w", local = TRUE)
  val = system2("install_name_tool", c("-id", filename, filename), stdout = NULL, stderr = NULL)
  if(val != 0) {
      e = simpleError("can't change install_name", call = NULL)
      stop(e)
  }
# print( getCurrentSOName(filename) )
  filename      
}


getCurrentSOName =
function(filename)
{
  txt = system(sprintf("otool -L %s", filename), intern = TRUE)
  ans = grep("^[[:space:]]+.*libnimble.dylib", txt, value = TRUE)
  trim(gsub("\\(compatibility.*", "", ans))
}

trim =
function (x) 
   gsub("(^[[:space:]]+|[[:space:]]+$)", "", x)
