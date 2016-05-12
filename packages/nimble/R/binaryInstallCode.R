
.onLoad =
function(...)
{
   if(!sameDir(system.file(package = "nimble"), AutoconfInfo$rpkg_install_dir)) {

       tryCatch(rebakeSOPath(libnimbleDLL()),
                 error = function(e) {
                            copyNimbleDylib()
                         })
   }
   
}


copyNimbleDylib =
function(tmpdir = tempdir())
{
    NimbleSessionInfo$soDir = tmpdir
    f = sprintf("%s/%s", tmpdir, basename(libnimbleDLL()))
    cat("Copying .so\n")       
    file.copy(libnimbleDLL(), tmpdir)
    # See if we have the write permissions to change the id in place.
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
  cmd = sprintf("install_name_tool -id  %s %s", filename, filename)

# print( getCurrentSOName(filename) )  
  system(cmd, intern = TRUE)
# print( getCurrentSOName(filename) )
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
